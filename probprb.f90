!  probprb.f90  11 July 2000
!
program probprb
!
!*******************************************************************************
!
!! PROBPRB calls sample problems for the PROB routines.
!
  write ( *, * ) ' '
  write ( *, * ) 'PROBPRB'
  write ( *, * ) '  Sample problems for the PROB routines.'

    call test201
    call test202

  call test001
  call test003
  call test004
  call test005
  call test006
  call test007
  call test008
  call test009
  call test010

  call test011
  call test013
  call test014
  call test015
  call test016
  call test017
  call test018
  call test019
  call test020

  call test021
  call test022

    call test220
    call test221

  call test023
  call test024
  call test026
  call test027
  call test028

    call test527

  call test030

  call test031
  call test033

    call test909

  call test034

    call test224
    call test225
    call test226

  call test035
  call test036
  call test037
  call test038
  call test040

  call test041
  call test042
  call test043
  call test044
  call test045
  call test046
  call test047
  call test049
  call test050

  call test052
  call test053
  call test054
  call test056
  call test057
  call test059
  call test060

  call test062
  call test063
  call test0645
  call test065
  call test066
  call test067
  call test069
  call test070

  call test072
  call test073
  call test074
  call test075

    call test620
    call test621

  call test077
  call test078
  call test079

    call test0795
    call test0796

    call test456
    call test457

  call test080

  call test081
  call test083
  call test085
  call test086
  call test087
  call test088
  call test089
  call test090

  call test092
  call test093
  call test095
  call test096
  call test098
  call test099
  call test100

  call test101
  call test102
  call test104
  call test105
  call test106
  call test107
  call test108

    call test520
    call test521

  call test109
  call test110

  call test112
  call test113
  call test115
  call test116
  call test118
  call test119
  call test120

  call test121
  call test122
  call test123
  call test124
  call test125
  call test126

    call test209
    call test211

  call test127
  call test128

    call test380
    call test381

    call test304
    call test305

    call test204
    call test205

  call test130

  call test131
  call test133
  call test134
  call test135
  call test137
  call test138
  call test139
  call test140

  call test142
  call test143
  call test144
  call test145
  call test146
  call test148
  call test149

  call test151
  call test152
  call test153
  call test154
  call test155

  write ( *, * ) ' '
  write ( *, * ) 'PROBPRB'
  write ( *, * ) '  Normal end of PROB tests.'

  stop
end
subroutine test201
!
!*******************************************************************************
!
!! TEST201 tests ANGLIT_CDF.
!! TEST201 tests ANGLIT_CDF_INV.
!! TEST201 tests ANGLIT_PDF.
!
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST201'
  write ( *, * ) '  For the Anglit PDF:'
  write ( *, * ) '  ANGLIT_CDF evaluates the CDF;'
  write ( *, * ) '  ANGLIT_CDF_INV inverts the CDF.'
  write ( *, * ) '  ANGLIT_PDF evaluates the PDF;'

  x = 0.50

  call anglit_pdf ( x, pdf )

  call anglit_cdf ( x, cdf )

  call anglit_cdf_inv ( cdf, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test202
!
!*******************************************************************************
!
!! TEST202 tests ANGLIT_MEAN;
!! TEST202 tests ANGLIT_SAMPLE;
!! TEST202 tests ANGLIT_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST202'
  write ( *, * ) '  For the Anglit PDF:'
  write ( *, * ) '  ANGLIT_MEAN computes the mean;'
  write ( *, * ) '  ANGLIT_SAMPLE samples;'
  write ( *, * ) '  ANGLIT_VARIANCE computes the variance.'

  call get_seed ( iseed )
  call anglit_mean ( mean )
  call anglit_variance ( variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call anglit_sample ( iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test001
!
!*******************************************************************************
!
!! TEST001 tests ARCSIN_CDF.
!! TEST001 tests ARCSIN_CDF_INV.
!! TEST001 tests ARCSIN_PDF.
!
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST001'
  write ( *, * ) '  For the Arcsin PDF:'
  write ( *, * ) '  ARCSIN_CDF evaluates the CDF;'
  write ( *, * ) '  ARCSIN_CDF_INV inverts the CDF.'
  write ( *, * ) '  ARCSIN_PDF evaluates the PDF;'

  x = 0.50

  call arcsin_pdf ( x, pdf )

  call arcsin_cdf ( x, cdf )

  call arcsin_cdf_inv ( cdf, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test003
!
!*******************************************************************************
!
!! TEST003 tests ARCSIN_MEAN;
!! TEST003 tests ARCSIN_SAMPLE;
!! TEST003 tests ARCSIN_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST003'
  write ( *, * ) '  For the Arcsin PDF:'
  write ( *, * ) '  ARCSIN_MEAN computes the mean;'
  write ( *, * ) '  ARCSIN_SAMPLE samples;'
  write ( *, * ) '  ARCSIN_VARIANCE computes the variance.'

  call get_seed ( iseed )
  call arcsin_mean ( mean )
  call arcsin_variance ( variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call arcsin_sample ( iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test004
!
!*******************************************************************************
!
!! TEST004 tests BENFORD_PDF.
!
  integer n
  real pdf
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST004'
  write ( *, * ) '  For the Benford PDF:'
  write ( *, * ) '  BENFORD_PDF evaluates the PDF.'

  write ( *, * ) ' '
  write ( *, * ) '  N    PDF(N)'
  write ( *, * ) ' '

  do n = 1, 19

    call benford_pdf ( n, pdf )
    write ( *, '(i6,g14.6)' ) n, pdf

  end do

  return
end
subroutine test005
!
!*******************************************************************************
!
!! TEST005 tests BERNOULLI_CDF.
!! TEST005 tests BERNOULLI_CDF_INV.
!
  real a
  real cdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST005'
  write ( *, * ) '  For the Bernoulli PDF,'
  write ( *, * ) '  BERNOULLI_CDF evaluates the CDF;'
  write ( *, * ) '  BERNOULLI_CDF_INV inverts the CDF.'

  cdf = 0.5
  a = 0.75

  call bernoulli_check ( a )

  call bernoulli_cdf_inv ( cdf, a, x )

  write ( *, * ) ' '
  write ( *, * ) '  CDF_INV argument CDF =        ', cdf
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  CDF_INV value X =             ', X
  write ( *, * ) ' '
  write ( *, * ) '  (Expected answer is 0)'

  call bernoulli_cdf ( x, a, cdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) ' '
  write ( *, * ) '  (Expected answer is 0.75)'

  return
end
subroutine test006
!
!*******************************************************************************
!
!! TEST006 tests BERNOULLI_MEAN;
!! TEST006 tests BERNOULLI_SAMPLE;
!! TEST006 tests BERNOULLI_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST006'
  write ( *, * ) '  For the Bernoulli PDF:'
  write ( *, * ) '  BERNOULLI_MEAN computes the mean;'
  write ( *, * ) '  BERNOULLI_SAMPLE samples;'
  write ( *, * ) '  BERNOULLI_VARIANCE computes the variance.'

  a = 0.75
  call bernoulli_check ( a )

  call get_seed ( iseed )
  call bernoulli_mean ( a, mean )
  call bernoulli_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call bernoulli_sample ( a, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test007
!
!*******************************************************************************
!
!! TEST007 tests BERNOULLI_PDF.
!! TEST007 tests BERNOULLI_CDF.
!
  real a
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST007'
  write ( *, * ) '  For the Bernoulli PDF:'
  write ( *, * ) '  BERNOULLI_PDF evaluates the PDF.'
  write ( *, * ) '  BERNOULLI_CDF evaluates the CDF.'

  x = 1

  a = 0.75

  call bernoulli_check ( a )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) ' '
  write ( *, * ) '     X      PDF(X)      CDF(X)'
  write ( *, * ) ' '

  do x = 0, 1
    call bernoulli_pdf ( x, a, pdf )
    call bernoulli_cdf ( x, a, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test008
!
!*******************************************************************************
!
!! TEST008 tests BETA;
!! TEST008 tests GAMMA.
!
  real a
  real b
  real beta
  real beta1
  real beta2
  real gamma
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST008'
  write ( *, * ) '  BETA evaluates the Beta function;'
  write ( *, * ) '  GAMMA evaluates the Gamma function.'

  a = 2.2
  b = 3.7

  beta1 = beta ( a, b )
  beta2 = gamma ( a ) * gamma ( b ) / gamma ( a + b )

  write ( *, * ) ' '
  write ( *, * ) '  Argument A =                   ', a
  write ( *, * ) '  Argument B =                   ', b
  write ( *, * ) '  Beta(A,B) =                    ', beta1
  write ( *, * ) '  (Expected value = 0.0454 )'
  write ( *, * ) ' '
  write ( *, * ) '  Gamma(A)*Gamma(B)/Gamma(A+B) = ', beta2

  return
end
subroutine test009
!
!*******************************************************************************
!
!! TEST009 tests BETA_CDF;
!! TEST009 tests BETA_CDF_INV.
!! TEST009 tests BETA_PDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST009'
  write ( *, * ) '  For the Beta PDF:'
  write ( *, * ) '  BETA_CDF evaluates the CDF;'
  write ( *, * ) '  BETA_CDF_INV inverts the CDF.'
  write ( *, * ) '  BETA_PDF evaluates the PDF;'

  x = 0.6

  a = 12.0
  b = 12.0

  call beta_check ( a, b )

  call beta_pdf ( x, a, b, pdf )

  call beta_cdf ( x, a, b, cdf )

  call beta_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =        ', x
  write ( *, * ) '  PDF parameter A =       ', a
  write ( *, * ) '  PDF parameter B =       ', b
  write ( *, * ) '  PDF value =             ', pdf
  write ( *, * ) '  CDF value =             ', cdf
  write ( *, * ) '  CDF_INV value X =       ', x2

  return
end
subroutine test010
!
!*******************************************************************************
!
!! TEST010 tests BETA_INC.
!
  real a
  real b
  real beta_inc
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST010'
  write ( *, * ) '  BETA_INC evaluates the incomplete Beta'
  write ( *, * ) '    function.'

  x = 0.61

  a = 2.2
  b = 3.7

  write ( *, * ) ' '
  write ( *, * ) '  X = ', x
  write ( *, * ) '  A = ', a
  write ( *, * ) '  B = ', b

  write ( *, * ) '  BETA_INC(X,A,B) = ', beta_inc ( x, a, b )
  write ( *, * ) '  (Expected value = 0.8822 )'

  return
end
subroutine test011
!
!*******************************************************************************
!
!! TEST011 tests BETA_MEAN;
!! TEST011 tests BETA_SAMPLE;
!! TEST011 tests BETA_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST011'
  write ( *, * ) '  For the Beta PDF:'
  write ( *, * ) '  BETA_MEAN computes the mean;'
  write ( *, * ) '  BETA_SAMPLE samples;'
  write ( *, * ) '  BETA_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0

  call beta_check ( a, b )

  call get_seed ( iseed )
  call beta_mean ( a, b, mean )
  call beta_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call beta_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test013
!
!*******************************************************************************
!
!! TEST013 tests BETA_BINOMIAL_CDF.
!! TEST013 tests BETA_BINOMIAL_CDF_INV.
!! TEST013 tests BETA_BINOMIAL_PDF.
!
  real a
  real b
  integer c
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST013'
  write ( *, * ) '  For the Beta Binomial PDF,'
  write ( *, * ) '  BETA_BINOMIAL_CDF evaluates the CDF;'
  write ( *, * ) '  BETA_BINOMIAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  BETA_BINOMIAL_PDF evaluates the PDF;'

  x = 2

  a = 2.0
  b = 3.0
  c = 4

  call beta_binomial_check ( a, b, c )

  call beta_binomial_pdf ( x, a, b, c, pdf )

  call beta_binomial_cdf ( x, a, b, c, cdf )

  call beta_binomial_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value =               ', x2

  return
end
subroutine test014
!
!*******************************************************************************
!
!! TEST014 tests BETA_BINOMIAL_MEAN;
!! TEST014 tests BETA_BINOMIAL_SAMPLE;
!! TEST014 tests BETA_BINOMIAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST014'
  write ( *, * ) '  For the Beta Binomial PDF:'
  write ( *, * ) '  BETA_BINOMIAL_MEAN computes the mean;'
  write ( *, * ) '  BETA_BINOMIAL_SAMPLE samples;'
  write ( *, * ) '  BETA_BINOMIAL_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0
  c = 4

  call beta_binomial_check ( a, b, c )

  call get_seed ( iseed )
  call beta_binomial_mean ( a, b, c, mean )
  call beta_binomial_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call beta_binomial_sample ( a, b, c, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test015
!
!*******************************************************************************
!
!! TEST015 tests BETA_BINOMIAL_CDF.
!! TEST015 tests BETA_BINOMIAL_PDF.
!
  real a
  real b
  integer c
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST015'
  write ( *, * ) '  For the Beta Binomial PDF:'
  write ( *, * ) '  BETA_BINOMIAL_PDF evaluates the PDF;'
  write ( *, * ) '  BETA_BINOMIAL_CDF evaluates the CDF.'

  a = 2.0
  b = 3.0
  c = 4

  call beta_binomial_check ( a, b, c )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A =  ', a
  write ( *, * ) '  PDF parameter B =  ', b
  write ( *, * ) '  PDF parameter C =  ', c

  write ( *, * ) ' '
  write ( *, * ) '     X     PDF(X)     CDF(X)'
  write ( *, * ) ' '

  do x = 0, c

    call beta_binomial_pdf ( x, a, b, c, pdf )
    call beta_binomial_cdf ( x, a, b, c, cdf )

    write ( *, '(i6,2g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test016
!
!*******************************************************************************
!
!! TEST016 tests BETA_PASCAL_CDF.
!! TEST016 tests BETA_PASCAL_CDF_INV.
!! TEST016 tests BETA_PASCAL_PDF.
!
  integer a
  real b
  real c
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST016'
  write ( *, * ) '  For the Beta Pascal PDF:'
  write ( *, * ) '  BETA_PASCAL_CDF evaluates the CDF;'
  write ( *, * ) '  BETA_PASCAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  BETA_PASCAL_PDF evaluates the PDF;'

  a = 5
  b = 3.0
  c = 4.0

  x = a + 2

  call beta_pascal_check ( a, b, c )

  call beta_pascal_pdf ( x, a, b, c, pdf )

  call beta_pascal_cdf ( x, a, b, c, cdf )

  call beta_pascal_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value =               ', x2

  return
end
subroutine test017
!
!*******************************************************************************
!
!! TEST017 tests BETA_PASCAL_MEAN;
!! TEST017 tests BETA_PASCAL_SAMPLE;
!! TEST017 tests BETA_PASCAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST017'
  write ( *, * ) '  For the Beta Pascal PDF:'
  write ( *, * ) '  BETA_PASCAL_MEAN computes the mean;'
  write ( *, * ) '  BETA_PASCAL_SAMPLE samples;'
  write ( *, * ) '  BETA_PASCAL_VARIANCE computes the variance.'

  a = 5
  b = 3.0
  c = 4.0

  call beta_pascal_check ( a, b, c )

  call get_seed ( iseed )
  call beta_pascal_mean ( a, b, c, mean )
  call beta_pascal_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  write ( *, * ) ' '
  write ( *, * ) 'TEST017 - DEBUG - '
  write ( *, * ) '  BETA_PASCAL_SAMPLE is still goofy!'
  return

  do i = 1, nsample
    call beta_pascal_sample ( a, b, c, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test018
!
!*******************************************************************************
!
!! TEST018 tests BETA_PASCAL_CDF.
!! TEST018 tests BETA_PASCAL_PDF.
!
  integer a
  real b
  real c
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST018'
  write ( *, * ) '  For the Beta Pascal PDF:'
  write ( *, * ) '  BETA_PASCAL_PDF evaluates the PDF;'
  write ( *, * ) '  BETA_PASCAL_CDF evaluates the CDF.'

  a = 5
  b = 3.0
  c = 4.0

  call beta_pascal_check ( a, b, c )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A =  ', a
  write ( *, * ) '  PDF parameter B =  ', b
  write ( *, * ) '  PDF parameter C =  ', c

  write ( *, * ) ' '
  write ( *, * ) '     X      PDF(X)      CDF(X)'
  write ( *, * ) ' '

  do x = a, a + 10

    call beta_pascal_pdf ( x, a, b, c, pdf )
    call beta_pascal_cdf ( x, a, b, c, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test019
!
!*******************************************************************************
!
!! TEST019 tests BINOMIAL_CDF;
!! TEST019 tests BINOMIAL_CDF_INV.
!! TEST019 tests BINOMIAL_PDF;
!
  integer a
  real b
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST019'
  write ( *, * ) '  For the Binomial PDF:'
  write ( *, * ) '  BINOMIAL_CDF evaluates the CDF;'
  write ( *, * ) '  BINOMIAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  BINOMIAL_PDF evaluates the PDF;'

  x = 3

  a = 5
  b = 0.95

  call binomial_check ( a, b )

  call binomial_pdf ( x, a, b, pdf )

  call binomial_cdf ( x, a, b, cdf )

  call binomial_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =        ', x
  write ( *, * ) '  PDF parameter A =       ', a
  write ( *, * ) '  PDF parameter B =       ', b
  write ( *, * ) '  PDF value =             ', pdf
  write ( *, * ) '  CDF value =             ', cdf
  write ( *, * ) '  CDF_INV value X =       ', x2

  return
end
subroutine test020
!
!*******************************************************************************
!
!! TEST020 tests BINOMIAL_COEF;
!! TEST020 tests BINOMIAL_COEF_LOG.
!
  integer cnk1
  real cnk2_log
  real cnk2
  integer k
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST020'
  write ( *, * ) '  BINOMIAL_COEF evaluates binomial coefficients.'
  write ( *, * ) '  BINOMIAL_COEF_LOG evaluates the logarithm.'

  write ( *, * ) ' '
  write ( *, * ) 'N, K, C(N,K)'
  write ( *, * ) ' '

  do n = 0, 4
    do k = 0, n
      call binomial_coef ( n, k, cnk1 )
      call binomial_coef_log ( n, k, cnk2_log )
      cnk2 = exp ( cnk2_log )
      write ( *, '(3i6,g14.6)' ) n, k, cnk1, cnk2
    end do
  end do

  return
end
subroutine test021
!
!*******************************************************************************
!
!! TEST021 tests BINOMIAL_MEAN;
!! TEST021 tests BINOMIAL_SAMPLE;
!! TEST021 tests BINOMIAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST021'
  write ( *, * ) '  For the Binomial PDF:'
  write ( *, * ) '  BINOMIAL_MEAN computes the mean;'
  write ( *, * ) '  BINOMIAL_SAMPLE samples;'
  write ( *, * ) '  BINOMIAL_VARIANCE computes the variance.'

  a = 5
  b = 0.30

  call binomial_check ( a, b )

  call get_seed ( iseed )
  call binomial_mean ( a, b, mean )
  call binomial_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call binomial_sample ( a, b, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test022
!
!*******************************************************************************
!
!! TEST022 tests BINOMIAL_CDF.
!! TEST022 tests BINOMIAL_PDF.
!
  integer a
  real b
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST022'
  write ( *, * ) '  For the Binomial PDF:'
  write ( *, * ) '  BINOMIAL_PDF evaluates the PDF.'
  write ( *, * ) '  BINOMIAL_CDF evaluates the CDF.'

  a = 5
  b = 0.95

  call binomial_check ( a, b )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) '  PDF parameter B = ', b
  write ( *, * ) ' '
  write ( *, * ) '  X     PDF(X)    CDF(X)'
  write ( *, * ) ' '

  do x = -1, a + 1

    call binomial_pdf ( x, a, b, pdf )
    call binomial_cdf ( x, a, b, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test220
!
!*******************************************************************************
!
!! TEST220 tests BRADFORD_CDF.
!! TEST220 tests BRADFORD_CDF_INV.
!! TEST220 tests BRADFORD_PDF.
!
  real a
  real b
  real c
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST220'
  write ( *, * ) '  For the Bradford PDF:'
  write ( *, * ) '  BRADFORD_CDF evaluates the CDF;'
  write ( *, * ) '  BRADFORD_CDF_INV inverts the CDF.'
  write ( *, * ) '  BRADFORD_PDF evaluates the PDF;'

  x = 1.25

  a = 1.0
  b = 2.0
  c = 3.0

  call bradford_check ( a, b, c )

  call bradford_pdf ( x, a, b, c, pdf )

  call bradford_cdf ( x, a, b, c, cdf )

  call bradford_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test221
!
!*******************************************************************************
!
!! TEST221 tests BRADFORD_MEAN;
!! TEST221 tests BRADFORD_SAMPLE;
!! TEST221 tests BRADFORD_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST221'
  write ( *, * ) '  For the Bradford PDF:'
  write ( *, * ) '  BRADFORD_MEAN computes the mean;'
  write ( *, * ) '  BRADFORD_SAMPLE samples;'
  write ( *, * ) '  BRADFORD_VARIANCE computes the variance.'

  a = 1.0
  b = 2.0
  c = 3.0

  call bradford_check ( a, b, c )

  call get_seed ( iseed )
  call bradford_mean ( a, b, c, mean )
  call bradford_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call bradford_sample ( a, b, c, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test023
!
!*******************************************************************************
!
!! TEST023 tests BURR_CDF.
!! TEST023 tests BURR_CDF_INV.
!! TEST023 tests BURR_PDF.
!
  real a
  real b
  real c
  real cdf
  real d
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST023'
  write ( *, * ) '  For the Burr PDF:'
  write ( *, * ) '  BURR_CDF evaluates the CDF;'
  write ( *, * ) '  BURR_CDF_INV inverts the CDF.'
  write ( *, * ) '  BURR_PDF evaluates the PDF;'

  x = 3.0

  a = 1.0
  b = 2.0
  c = 3.0
  d = 2.0

  call burr_check ( a, b, c, d )

  call burr_pdf ( x, a, b, c, d, pdf )

  call burr_cdf ( x, a, b, c, d, cdf )

  call burr_cdf_inv ( cdf, a, b, c, d, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF parameter D =             ', d
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test024
!
!*******************************************************************************
!
!! TEST024 tests BURR_MEAN;
!! TEST024 tests BURR_VARIANCE;
!! TEST024 tests BURR_SAMPLE;
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  real d
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST024'
  write ( *, * ) '  For the Burr PDF:'
  write ( *, * ) '  BURR_MEAN computes the mean;'
  write ( *, * ) '  BURR_VARIANCE computes the variance;'
  write ( *, * ) '  BURR_SAMPLE samples;'

  a = 1.0
  b = 2.0
  c = 3.0
  d = 2.0

  call burr_check ( a, b, c, d )

  call get_seed ( iseed )

  call burr_mean ( a, b, c, d, mean )
  call burr_variance ( a, b, c, d, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF parameter C =         ', c
  write ( *, * ) '  PDF parameter D =         ', d
  write ( *, * ) '  PDF mean =                ', mean
  write ( *, * ) '  PDF variance =            ', variance

  do i = 1, nsample
    call burr_sample ( a, b, c, d, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test026
!
!*******************************************************************************
!
!! TEST026 tests CAUCHY_CDF.
!! TEST026 tests CAUCHY_CDF_INV.
!! TEST026 tests CAUCHY_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST026'
  write ( *, * ) '  For the Cauchy PDF:'
  write ( *, * ) '  CAUCHY_CDF evaluates the CDF;'
  write ( *, * ) '  CAUCHY_CDF_INV inverts the CDF.'
  write ( *, * ) '  CAUCHY_PDF evaluates the PDF;'

  x = 0.75

  a = 2.0
  b = 3.0

  call cauchy_check ( a, b )

  call cauchy_pdf ( x, a, b, pdf )

  call cauchy_cdf ( x, a, b, cdf )

  call cauchy_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test027
!
!*******************************************************************************
!
!! TEST027 tests CAUCHY_MEAN;
!! TEST027 tests CAUCHY_SAMPLE;
!! TEST027 tests CAUCHY_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST027'
  write ( *, * ) '  For the Cauchy PDF:'
  write ( *, * ) '  CAUCHY_MEAN computes the mean;'
  write ( *, * ) '  CAUCHY_VARIANCE computes the variance;'
  write ( *, * ) '  CAUCHY_SAMPLE samples.'

  a = 2.0
  b = 3.0

  call cauchy_check ( a, b )

  call get_seed ( iseed )
  call cauchy_mean ( a, b, mean )
  call cauchy_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF mean =                    ', variance
  
  do i = 1, nsample
    call cauchy_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test028
!
!*******************************************************************************
!
!! TEST028 tests CHI_CDF.
!! TEST028 tests CHI_CDF_INV.
!! TEST028 tests CHI_PDF.
!
  real a
  real b
  real c
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST028'
  write ( *, * ) '  For the Chi PDF:'
  write ( *, * ) '  CHI_CDF evaluates the CDF.'
  write ( *, * ) '  CHI_CDF_INV inverts the CDF.'
  write ( *, * ) '  CHI_PDF evaluates the PDF.'

  x = 2.0

  a = 1.0
  b = 2.0
  c = 3.0

  call chi_check ( a, b, c )

  call chi_pdf ( x, a, b, c, pdf )

  call chi_cdf ( x, a, b, c, cdf )

  call chi_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value =               ', x2

  return
end
subroutine test527
!
!*******************************************************************************
!
!! TEST527 tests CHI_MEAN;
!! TEST527 tests CHI_SAMPLE;
!! TEST527 tests CHI_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST527'
  write ( *, * ) '  For the Chi PDF:'
  write ( *, * ) '  CHI_MEAN computes the mean;'
  write ( *, * ) '  CHI_VARIANCE computes the variance;'
  write ( *, * ) '  CHI_SAMPLE samples.'

  a = 1.0
  b = 2.0
  c = 3.0

  call chi_check ( a, b, c )

  call get_seed ( iseed )
  call chi_mean ( a, b, c, mean )
  call chi_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF mean =                    ', variance
  
  do i = 1, nsample
    call chi_sample ( a, b, c, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test030
!
!*******************************************************************************
!
!! TEST030 tests CHISQUARE_CENTRAL_CDF.
!! TEST030 tests CHISQUARE_CENTRAL_CDF_INV.
!! TEST030 tests CHISQUARE_CENTRAL_PDF.
!
  real a
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST030'
  write ( *, * ) '  For the central chi square PDF:'
  write ( *, * ) '  CHISQUARE_CENTRAL_CDF evaluates the CDF;'
  write ( *, * ) '  CHISQUARE_CENTRAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  CHISQUARE_CENTRAL_PDF evaluates the PDF;'

  x = 6.0

  a = 4.0

  call chisquare_central_check ( a )

  call chisquare_central_pdf ( x, a, pdf )

  call chisquare_central_cdf ( x, a, cdf )

  call chisquare_central_cdf_inv ( cdf, a, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test031
!
!*******************************************************************************
!
!! TEST031 tests CHISQUARE_CENTRAL_MEAN;
!! TEST031 tests CHISQUARE_CENTRAL_SAMPLE;
!! TEST031 tests CHISQUARE_CENTRAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST031'
  write ( *, * ) '  For the central chi square PDF:'
  write ( *, * ) '  CHISQUARE_CENTRAL_MEAN computes the mean;'
  write ( *, * ) '  CHISQUARE_CENTRAL_SAMPLE samples;'
  write ( *, * ) '  CHISQUARE_CENTRAL_VARIANCE computes the variance.'

  a = 10.0

  call chisquare_central_check ( a )
  call get_seed ( iseed )
  call chisquare_central_mean ( a, mean )
  call chisquare_central_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call chisquare_central_sample ( a, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test033
!
!*******************************************************************************
!
!! TEST033 tests CHISQUARE_NONCENTRAL_MEAN;
!! TEST033 tests CHISQUARE_NONCENTRAL_SAMPLE;
!! TEST033 tests CHISQUARE_NONCENTRAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST033'
  write ( *, * ) '  For the noncentral chi square PDF:'
  write ( *, * ) '  CHISQUARE_NONCENTRAL_SAMPLE samples.'

  a = 3.0
  b = 2.0

  call chisquare_noncentral_check ( a, b )

  call get_seed ( iseed )
  call chisquare_noncentral_mean ( a, b, mean )
  call chisquare_noncentral_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call chisquare_noncentral_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test909
!
!*******************************************************************************
!
!! TEST909 tests CIRCLE_SAMPLE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  integer j
  real mean(2)
  real variance(2)
  real x_table(nsample,2)
  real x1
  real x2
  real xmax(2)
  real xmin(2)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST909'
  write ( *, * ) '  CIRCLE_SAMPLE samples points in a circle.'

  call get_seed ( iseed )

  a = 10.0
  b = 4.0
  c = 3.0

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  X coordinate of center is A = ', a
  write ( *, * ) '  Y coordinate of center is B = ', b
  write ( *, * ) '  Radius is C =                 ', c

  do i = 1, nsample
    call circle_sample ( a, b, c, iseed, x1, x2 )
    x_table(i,1) = x1
    x_table(i,2) = x2
  end do

  do j = 1, 2
    call rvec_mean ( nsample, x_table(1,j), mean(j) )
    call rvec_variance ( nsample, x_table(1,j), variance(j) )
    call rvec_max ( nsample, x_table(1,j), imax, xmax(j) )
    call rvec_min ( nsample, x_table(1,j), imin, xmin(j) )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean(1), mean(2)
  write ( *, * ) '  Sample variance = ', variance(1), variance(2)
  write ( *, * ) '  Sample maximum =  ', xmax(1), xmax(2)
  write ( *, * ) '  Sample minimum =  ', xmin(1), xmin(2)

  return
end
subroutine test034
!
!*******************************************************************************
!
!! TEST034 tests CIRCULAR_NORMAL_01_MEAN;
!! TEST034 tests CIRCULAR_NORMAL_01_SAMPLE;
!! TEST034 tests CIRCULAR_NORMAL_01_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer iseed
  integer imax
  integer imin
  integer j
  real mean(2)
  real variance(2)
  real x(2)
  real x_table(nsample,2)
  real xmax(2)
  real xmin(2)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST034'
  write ( *, * ) '  For the Circular Normal 01 PDF:'
  write ( *, * ) '  CIRCULAR_NORMAL_01_MEAN computes the mean;'
  write ( *, * ) '  CIRCULAR_NORMAL_01_SAMPLE samples;'
  write ( *, * ) '  CIRCULAR_NORMAL_01_VARIANCE computes variance.'

  call get_seed ( iseed )
  call circular_normal_01_mean ( mean )
  call circular_normal_01_variance ( variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF means =               ', mean(1), mean(2)
  write ( *, * ) '  PDF variances =           ', variance(1), variance(2)
  
  do i = 1, nsample
    call circular_normal_01_sample ( iseed, x )
    x_table(i,1) = x(1)
    x_table(i,2) = x(2)
  end do

  do j = 1, 2
    call rvec_mean ( nsample, x_table(1,j), mean(j) )
    call rvec_variance ( nsample, x_table(1,j), variance(j) )
    call rvec_max ( nsample, x_table(1,j), imax, xmax(j) )
    call rvec_min ( nsample, x_table(1,j), imin, xmin(j) )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean(1), mean(2)
  write ( *, * ) '  Sample variance = ', variance(1), variance(2)
  write ( *, * ) '  Sample maximum =  ', xmax(1), xmax(2)
  write ( *, * ) '  Sample minimum =  ', xmin(1), xmin(2)

  return
end
subroutine test224
!
!*******************************************************************************
!
!! TEST224 tests COSINE_CDF.
!! TEST224 tests COSINE_CDF_INV.
!! TEST224 tests COSINE_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST224'
  write ( *, * ) '  For the Cosine PDF:'
  write ( *, * ) '  COSINE_CDF evaluates the CDF.'
  write ( *, * ) '  COSINE_CDF_INV inverts the CDF.'
  write ( *, * ) '  COSINE_PDF evaluates the PDF.'

  x = 1.0

  a = 2.0
  b = 1.0

  call cosine_check ( a, b )

  call cosine_cdf ( x, a, b, cdf )

  call cosine_pdf ( x, a, b, pdf )

  call cosine_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value             ', x2

  return
end
subroutine test225
!
!*******************************************************************************
!
!! TEST225 tests COSINE_MEAN;
!! TEST225 tests COSINE_SAMPLE;
!! TEST225 tests COSINE_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST225'
  write ( *, * ) '  For the Cosine PDF:'
  write ( *, * ) '  COSINE_MEAN computes the mean;'
  write ( *, * ) '  COSINE_SAMPLE samples;'
  write ( *, * ) '  COSINE_VARIANCE computes the variance.'

  a = 2.0
  b = 1.0

  call cosine_check ( a, b )

  call get_seed ( iseed )
  call cosine_mean ( a, b, mean )
  call cosine_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call cosine_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test226
!
!*******************************************************************************
!
!! TEST226 tests COUPON_SIMULATE.
!
  integer, parameter :: n_trial = 10
  integer, parameter :: max_type = 25
!
  real average
  integer coupon(max_type)
  real expect
  integer i
  integer iseed
  integer n_coupon
  integer n_type
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST226'
  write ( *, * ) '  COUPON_SIMULATE simulates the coupon collector''s problem.'
  write ( *, * ) ' '

  call get_seed ( iseed )
  write ( *, * ) '  Using random seed ISEED = ', iseed

  do n_type = 5, 25, 5

    write ( *, * ) ' '
    write ( *, * ) '  Number of coupon types is ', n_type
    expect = real ( n_type ) * log ( real ( n_type ) )
    write ( *, * ) '  Expected wait is about ', expect
    write ( *, * ) ' '

    average = 0.0
    do i = 1, n_trial
      call coupon_simulate ( n_type, iseed, coupon, n_coupon )
      write ( *, '(2i5)' ) i, n_coupon
      average = average + real ( n_coupon )
    end do

    average = average / real ( n_trial )
    write ( *, * ) ' '
    write ( *, * ) '  Average wait was ', average

  end do

  return
end
subroutine test035
!
!*******************************************************************************
!
!! TEST035 tests DERANGED_CDF;
!! TEST035 tests DERANGED_CDF_INV.
!! TEST035 tests DERANGED_PDF;
!
  integer a
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST035'
  write ( *, * ) '  For the Deranged PDF:'
  write ( *, * ) '  DERANGED_CDF evaluates the CDF;'
  write ( *, * ) '  DERANGED_CDF_INV inverts the CDF.'
  write ( *, * ) '  DERANGED_PDF evaluates the PDF;'

  x = 3

  a = 7

  call deranged_check ( a )

  call deranged_pdf ( x, a, pdf )

  call deranged_cdf ( x, a, cdf )

  call deranged_cdf_inv ( cdf, a, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =        ', x
  write ( *, * ) '  PDF parameter A =       ', a
  write ( *, * ) '  PDF value =             ', pdf
  write ( *, * ) '  CDF value =             ', cdf
  write ( *, * ) '  CDF_INV value X =       ', x2

  return
end
subroutine test036
!
!*******************************************************************************
!
!! TEST036 tests DERANGED_CDF.
!! TEST036 tests DERANGED_PDF.
!
  integer a
  real cdf
  real pdf
  integer x
!
  a = 7

  write ( *, * ) ' '
  write ( *, * ) 'TEST036'
  write ( *, * ) '  For the Deranged PDF:'
  write ( *, * ) '  DERANGED_PDF evaluates the PDF.'
  write ( *, * ) '  DERANGED_CDF evaluates the CDF.'
  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) ' '
  write ( *, * ) '  X    PDF(X)      CDF(X)'
  write ( *, * ) ' '

  call deranged_check ( a )

  do x = 0, a
    call deranged_pdf ( x, a, pdf )
    call deranged_cdf ( x, a, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test037
!
!*******************************************************************************
!
!! TEST037 tests DERANGED_MEAN.
!! TEST037 tests DERANGED_VARIANCE.
!! TEST037 tests DERANGED_SAMPLE.
!
  integer, parameter :: nsample = 1000
!
  integer a
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST037'
  write ( *, * ) '  For the Deranged PDF:'
  write ( *, * ) '  DERANGED_MEAN computes the mean.'
  write ( *, * ) '  DERANGED_VARIANCE computes the variance.'
  write ( *, * ) '  DERANGED_SAMPLE samples.'

  call get_seed ( iseed )

  a = 7

  call deranged_check ( a )
  call deranged_mean ( a, mean )
  call deranged_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call deranged_sample ( a, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test038
!
!*******************************************************************************
!
!! TEST038 tests DIPOLE_CDF.
!! TEST038 tests DIPOLE_CDF_INV.
!! TEST038 tests DIPOLE_PDF.
!
  integer, parameter :: ntest = 3
!
  real a
  real atest(ntest)
  real b
  real btest(ntest)
  real cdf
  integer itest
  real pi
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST038'
  write ( *, * ) '  For the Dipole PDF:'
  write ( *, * ) '  DIPOLE_CDF evaluates the CDF.'
  write ( *, * ) '  DIPOLE_CDF_INV inverts the CDF.'
  write ( *, * ) '  DIPOLE_PDF evaluates the PDF.'

  atest(1) = 0.0
  btest(1) = 1.0
  atest(2) = pi() / 4.0
  btest(2) = 0.5
  atest(3) = pi() / 2.0
  btest(3) = 0.0

  do itest = 1, ntest

    x = 0.6

    a = atest(itest)
    b = btest(itest)

    call dipole_check ( a, b )

    call dipole_pdf ( x, a, b, pdf )

    call dipole_cdf ( x, a, b, cdf )

    call dipole_cdf_inv ( cdf, a, b, x2 )

    write ( *, * ) ' '
    write ( *, * ) '  PDF argument X =        ', x
    write ( *, * ) '  PDF parameter A =       ', a
    write ( *, * ) '  PDF parameter B =       ', b
    write ( *, * ) '  PDF value =             ', pdf
    write ( *, * ) '  CDF value =             ', cdf
    write ( *, * ) '  CDF_INV value =         ', x2

  end do

  return
end
subroutine test040
!
!*******************************************************************************
!
!! TEST040 tests DIPOLE_SAMPLE.
!
  integer, parameter :: nsample = 1000
  integer, parameter :: ntest = 3
!
  real a
  real atest(ntest)
  real b
  real btest(ntest)
  integer i
  integer iseed
  integer itest
  integer imax
  integer imin
  real mean
  real pi
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST040'
  write ( *, * ) '  For the Dipole PDF:'
  write ( *, * ) '  DIPOLE_SAMPLE samples.'

  atest(1) = 0.0
  btest(1) = 1.0
  atest(2) = pi () / 4.0
  btest(2) = 0.5
  atest(3) = pi() / 2.0
  btest(3) = 0.0

  do itest = 1, ntest

    a = atest(itest)
    b = btest(itest)

    call dipole_check ( a, b )

    call get_seed ( iseed )

    write ( *, * ) ' '
    write ( *, * ) '  Initial random number seed is ', iseed
    write ( *, * ) '  PDF parameter A =             ', a
    write ( *, * ) '  PDF parameter B =             ', b
  
    do i = 1, nsample
      call dipole_sample ( a, b, iseed, x(i) )
    end do

    call rvec_mean ( nsample, x, mean )
    call rvec_variance ( nsample, x, variance )
    call rvec_max ( nsample, x, imax, xmax )
    call rvec_min ( nsample, x, imin, xmin )

    write ( *, * ) ' '
    write ( *, * ) '  Sample size =     ', nsample
    write ( *, * ) '  Sample mean =     ', mean
    write ( *, * ) '  Sample variance = ', variance
    write ( *, * ) '  Sample maximum =  ', xmax
    write ( *, * ) '  Sample minimum =  ', xmin

  end do

  return
end
subroutine test041
!
!*******************************************************************************
!
!! TEST041 tests DIRICHLET_MEAN;
!! TEST041 tests DIRICHLET_SAMPLE.
!! TEST041 tests DIRICHLET_VARIANCE.
!
  integer, parameter :: n = 3
  integer, parameter :: nsample = 1000
!
  real a(n)
  integer i
  integer iseed
  integer imax(n)
  integer imin(n)
  integer j
  real mean(n)
  real m2(n,n)
  real variance(n)
  real x(n,nsample)
  real xmax(n)
  real xmin(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST041'
  write ( *, * ) '  For the Dirichlet PDF:'
  write ( *, * ) '  DIRICHLET_SAMPLE samples;'
  write ( *, * ) '  DIRICHLET_MEAN computes the mean;'
  write ( *, * ) '  DIRICHLET_VARIANCE computes the variance.'

  a(1) = 0.250
  a(2) = 0.500
  a(3) = 1.250

  call dirichlet_check ( n, a )

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is   ', iseed
  write ( *, * ) '  Number of components N =        ', n
  write ( *, * ) '  PDF parameters A(I), I = 1 to N:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  call dirichlet_mean ( n, a, mean )

  call dirichlet_variance ( n, a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  PDF mean, variance:'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(i6,2g14.6)' ) i, mean(i), variance(i)
  end do

  call dirichlet_moment2 ( n, a, m2 )
  
  write ( *, * ) ' '
  write ( *, * ) '  Second moments:'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3g14.6)' ) m2(i,1:n)
  end do

  do i = 1, nsample
    call dirichlet_sample ( n, a, iseed, x(1,i) )
  end do

  call rrow_max ( n, n, nsample, x, imax, xmax )
  call rrow_min ( n, n, nsample, x, imin, xmin )
  call rrow_mean ( n, n, nsample, x, mean )
  call rrow_variance ( n, n, nsample, x, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size = ', nsample
  write ( *, * ) '  Observed Min, Max, Mean, Variance:'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(i6,4g14.6)' ) i, xmin(i), xmax(i), mean(i), variance(i)
  end do

  return
end
subroutine test042
!
!*******************************************************************************
!
!! TEST042 tests DIRICHLET_PDF.
!
  integer, parameter :: n = 3
!
  real a(n)
  integer i
  real pdf
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST042'
  write ( *, * ) '  For the Dirichlet PDF:'
  write ( *, * ) '  DIRICHLET_PDF evaluates the PDF.'

  x(1) = 0.500
  x(2) = 0.125
  x(3) = 0.375

  a(1) = 0.250
  a(2) = 0.500
  a(3) = 1.250

  call dirichlet_check ( n, a )

  write ( *, * ) '  Number of components N =        ', n
  write ( *, * ) '  PDF arguments X(I), I = 1 to N:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, x(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameters A(I), I = 1 to N:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, a(i)
  end do

  call dirichlet_pdf ( x, n, a, pdf )

  write ( *, * ) 'PDF value =           ', pdf

  return
end
subroutine test043
!
!*******************************************************************************
!
!! TEST043 tests DIRICHLET_MIX_MEAN;
!! TEST043 tests DIRICHLET_MIX_SAMPLE.
!
  integer, parameter :: comp_num = 2
  integer, parameter :: elem_num = 3
  integer, parameter :: sample_num = 1000

  integer, parameter :: elem_max = elem_num
!
  real a(elem_max,comp_num)
  integer comp
  integer comp_i
  real comp_weight(comp_num)
  integer elem_i
  integer iseed
  integer imax(elem_num)
  integer imin(elem_num)
  real mean(elem_num)
  integer sample_i
  real variance(elem_num)
  real x(elem_num,sample_num)
  real xmax(elem_num)
  real xmin(elem_num)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST043'
  write ( *, * ) '  For the Dirichlet Mixture PDF:'
  write ( *, * ) '  DIRICHLET_MIX_SAMPLE samples;'
  write ( *, * ) '  DIRICHLET_MIX_MEAN computes the mean;'

  a(1,1) = 0.250
  a(2,1) = 0.500
  a(3,1) = 1.250

  a(1,2) = 2.000
  a(2,2) = 0.000
  a(3,2) = 2.000

  comp_weight(1) = 1.0
  comp_weight(2) = 2.0

  call dirichlet_mix_check ( comp_num, elem_max, elem_num, a, comp_weight )

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is   ', iseed
  write ( *, * ) '  Number of elements ELEM_NUM =   ', elem_num
  write ( *, * ) '  Number of components COMP_NUM = ', comp_num
  write ( *, * ) '  PDF parameters A(ELEM,COMP):'
  write ( *, * ) ' '
  do elem_i = 1, elem_num
    write ( *, '(2g14.6)' ) ( a(elem_i,comp_i), comp_i = 1, comp_num )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Component weights:'
  write ( *, * ) ' '
  do comp_i = 1, comp_num
    write ( *, '(i6,g14.6)' ) comp_i, comp_weight(comp_i)
  end do

  call dirichlet_mix_mean ( comp_num, elem_max, elem_num, a, comp_weight, mean )

  write ( *, * ) ' '
  write ( *, * ) '  PDF mean:'
  write ( *, * ) ' '

  do elem_i = 1, elem_num
    write ( *, '(i6,2g14.6)' ) elem_i, mean(elem_i)
  end do

  do sample_i = 1, sample_num
    call dirichlet_mix_sample ( comp_num, elem_max, elem_num, a, &
      comp_weight, iseed, comp, x(1,sample_i) )
  end do

  call rrow_max ( elem_num, elem_num, sample_num, x, imax, xmax )

  call rrow_min ( elem_num, elem_num, sample_num, x, imin, xmin )

  call rrow_mean ( elem_num, elem_num, sample_num, x, mean )

  call rrow_variance ( elem_num, elem_num, sample_num, x, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size = ', sample_num
  write ( *, * ) '  Observed Min, Max, Mean, Variance:'
  write ( *, * ) ' '

  do elem_i = 1, elem_num
    write ( *, '(i6,4g14.6)' ) elem_i, xmin(elem_i), xmax(elem_i), &
      mean(elem_i), variance(elem_i)
  end do

  return
end
subroutine test044
!
!*******************************************************************************
!
!! TEST044 tests DIRICHLET_MIX_PDF.
!
  integer, parameter :: comp_num = 2
  integer, parameter :: elem_num = 3

  integer, parameter :: elem_max = elem_num
!
  real a(elem_max,comp_num)
  integer comp_i
  real comp_weight(comp_num)
  integer elem_i
  real pdf
  real x(elem_num)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST044'
  write ( *, * ) '  For the Dirichlet mixture PDF:'
  write ( *, * ) '  DIRICHLET_MIX_PDF evaluates the PDF.'

  x(1) = 0.500
  x(2) = 0.125
  x(3) = 0.375

  a(1,1) = 0.250
  a(2,1) = 0.500
  a(3,1) = 1.250

  a(1,2) = 2.000
  a(2,2) = 0.000
  a(3,2) = 2.000

  comp_weight(1) = 1.0
  comp_weight(2) = 2.0

  call dirichlet_mix_check ( comp_num, elem_max, elem_num, a, comp_weight )

  write ( *, * ) '  Number of elements ELEM_NUM =   ', elem_num
  write ( *, * ) '  Number of components COMP_NUM = ', comp_num
  write ( *, * ) '  PDF parameters A(ELEM,COMP):'
  write ( *, * ) ' '
  do elem_i = 1, elem_num
    write ( *, '(2g14.6)' ) a(elem_i,1:comp_num)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Component weights:'
  write ( *, * ) ' '
  do comp_i = 1, comp_num
    write ( *, '(i6,g14.6)' ) comp_i, comp_weight(comp_i)
  end do

  call dirichlet_mix_pdf ( x, comp_num, elem_max, elem_num, a, comp_weight, &
    pdf )

  write ( *, * ) 'PDF value =           ', pdf

  return
end
subroutine test045
!
!*******************************************************************************
!
!! TEST045 tests BETA_PDF.
!! TEST045 tests DIRICHLET_PDF.
!
  integer, parameter :: n = 2
!
  real a
  real aval
  real avec(n)
  real b
  real bval
  integer i
  real pdf
  real x
  real xval
  real xvec(n)
!
  xval = 0.25
  aval = 2.5
  bval = 3.5

  write ( *, * ) ' '
  write ( *, * ) 'TEST045'
  write ( *, * ) '  BETA_PDF evaluates the Beta PDF.'
  write ( *, * ) '  DIRICHLET_PDF evaluates the Dirichlet PDF.'
  write ( *, * ) ' '
  write ( *, * ) '  For N = 2, Dirichlet = Beta.'

  xvec(1) = xval
  xvec(2) = 1.0 - xval

  avec(1) = aval
  avec(2) = bval

  call dirichlet_check ( n, avec )

  write ( *, * ) '  Number of components N =        ', n
  write ( *, * ) '  PDF arguments X(I), I = 1 to N:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, xvec(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameters A(I), I = 1 to N:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, avec(i)
  end do

  call dirichlet_pdf ( xvec, n, avec, pdf )

  write ( *, * ) 'Dirichlet PDF value =  ', pdf

  x = xval

  a = aval
  b = bval

  call beta_pdf ( x, a, b, pdf )

  write ( *, * ) 'Beta PDF value =       ', pdf

  return
end
subroutine test046
!
!*******************************************************************************
!
!! TEST046 tests DISCRETE_CDF.
!! TEST046 tests DISCRETE_CDF_INV.
!! TEST046 tests DISCRETE_PDF.
!
  integer, parameter :: a = 6
!
  real b(a)
  real cdf
  integer j
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST046'
  write ( *, * ) '  For the Discrete PDF:'
  write ( *, * ) '  DISCRETE_CDF evaluates the CDF;'
  write ( *, * ) '  DISCRETE_CDF_INV inverts the CDF.'
  write ( *, * ) '  DISCRETE_PDF evaluates the PDF;'

  b(1) = 1.0
  b(2) = 2.0
  b(3) = 6.0
  b(4) = 2.0
  b(5) = 4.0
  b(6) = 1.0

  call discrete_check ( a, b )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) '  PDF parameters B:'
  do j = 1, a
    write ( *, '(i6,g14.6)' ) j, b(j)
  end do

  write ( *, * ) ' '
  write ( *, * ) '    X    PDF    CDF  CDF_INV(CDF)'
  write ( *, * ) ' '

  do x = 0, 7

    call discrete_pdf ( x, a, b, pdf )

    call discrete_cdf ( x, a, b, cdf )

    call discrete_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(i6,2g14.6,i6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test047
!
!*******************************************************************************
!
!! TEST047 tests DISCRETE_MEAN;
!! TEST047 tests DISCRETE_SAMPLE;
!! TEST047 tests DISCRETE_VARIANCE.
!
  integer, parameter :: a = 6
  integer, parameter :: nsample = 1000
!
  real b(a)
  integer i
  integer iseed
  integer imax
  integer imin
  integer j
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST047'
  write ( *, * ) '  For the Discrete PDF:'
  write ( *, * ) '  DISCRETE_MEAN computes the mean;'
  write ( *, * ) '  DISCRETE_SAMPLE samples;'
  write ( *, * ) '  DISCRETE_VARIANCE computes the variance.'

  b(1) = 1.0
  b(2) = 2.0
  b(3) = 6.0
  b(4) = 2.0
  b(5) = 4.0
  b(6) = 1.0

  call discrete_check ( a, b )
  call get_seed ( iseed )
  call discrete_mean ( a, b, mean )
  call discrete_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameters B:'
  do j = 1, a
    write ( *, '(i6,g14.6)' ) j, b(j)
  end do
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call discrete_sample ( a, b, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test049
!
!*******************************************************************************
!
!! TEST049 tests ERLANG_CDF.
!! TEST049 tests ERLANG_CDF_INV.
!! TEST049 tests ERLANG_PDF.
!
  real a
  real b
  integer c
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST049'
  write ( *, * ) '  For the Erlang PDF:'
  write ( *, * ) '  ERLANG_CDF evaluates the CDF.'
  write ( *, * ) '  ERLANG_CDF_INV inverts the CDF.'
  write ( *, * ) '  ERLANG_PDF evaluates the PDF.'

  x = 4.0

  a = 1.0
  b = 2.0
  c = 3

  call erlang_check ( a, b, c )

  call erlang_pdf ( x, a, b, c, pdf )

  call erlang_cdf ( x, a, b, c, cdf )

  call erlang_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF parameter C =         ', c
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV                   ', x2

  return
end
subroutine test050
!
!*******************************************************************************
!
!! TEST050 tests ERLANG_MEAN;
!! TEST050 tests ERLANG_SAMPLE;
!! TEST050 tests ERLANG_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST050'
  write ( *, * ) '  For the Erlang PDF:'
  write ( *, * ) '  ERLANG_MEAN computes the mean;'
  write ( *, * ) '  ERLANG_SAMPLE samples;'
  write ( *, * ) '  ERLANG_VARIANCE computes the variance.'

  a = 1.0
  b = 2.0
  c = 3

  call erlang_check ( a, b, c )

  call get_seed ( iseed )
  call erlang_mean ( a, b, c, mean )
  call erlang_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call erlang_sample ( a, b, c, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test052
!
!*******************************************************************************
!
!! TEST052 tests ERROR_FUNCTION.
!
  real cdf
  real erfx
  real error_function
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST052'
  write ( *, * ) '  ERROR_FUNCTION evaluates ERF(X).'

  x = 1.0

  erfx = error_function ( x )

  write ( *, * ) ' '
  write ( *, * ) '  ERF argument X =          ', x
  write ( *, * ) '  ERF value                 ', erfx
  write ( *, * ) ' '
  write ( *, * ) '  (Expected answer is 0.843)'
  write ( *, * ) ' '
  write ( *, * ) '  Test:'
  write ( *, * ) '    0.5 * ( ERF(X/SQRT(2)) + 1 ) = Normal_CDF(X)'
  write ( *, * ) ' '

  x = 1.0
  x2 = x / sqrt ( 2.0 )
  erfx = error_function ( x2 )

  call normal_01_cdf ( x, cdf )

  write ( *, * ) '0.5 * ( ERF(X/SQRT(2)) + 1 ) = ', 0.5 * ( erfx + 1.0 )
  write ( *, * ) 'Normal_CDF(X) = ', cdf

  return
end
subroutine test053
!
!*******************************************************************************
!
!! TEST053 tests EXPONENTIAL_01_CDF;
!! TEST053 tests EXPONENTIAL_01_CDF_INV.
!! TEST053 tests EXPONENTIAL_01_PDF;
!
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST053'
  write ( *, * ) '  For the Exponential 01 PDF:'
  write ( *, * ) '  EXPONENTIAL_01_CDF evaluates the CDF.'
  write ( *, * ) '  EXPONENTIAL_01_CDF_INV inverts the CDF.'
  write ( *, * ) '  EXPONENTIAL_01_PDF evaluates the PDF.'

  x = 0.5

  call exponential_01_pdf ( x, pdf )

  call exponential_01_cdf ( x, cdf )

  call exponential_01_cdf_inv ( cdf, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value X           ', x2

  return
end
subroutine test054
!
!*******************************************************************************
!
!! TEST054 tests EXPONENTIAL_01_MEAN;
!! TEST054 tests EXPONENTIAL_01_SAMPLE;
!! TEST054 tests EXPONENTIAL_01_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST054'
  write ( *, * ) '  For the Exponential 01_PDF:'
  write ( *, * ) '  EXPONENTIAL_01_MEAN computes the mean;'
  write ( *, * ) '  EXPONENTIAL_01_SAMPLE samples;'
  write ( *, * ) '  EXPONENTIAL_01_VARIANCE computes the variance.'

  call get_seed ( iseed )
  call exponential_01_mean ( mean )
  call exponential_01_variance ( variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call exponential_01_sample ( iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test056
!
!*******************************************************************************
!
!! TEST056 tests EXPONENTIAL_CDF;
!! TEST056 tests EXPONENTIAL_CDF_INV.
!! TEST056 tests EXPONENTIAL_PDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST056'
  write ( *, * ) '  For the Exponential CDF:'
  write ( *, * ) '  EXPONENTIAL_CDF evaluates the CDF.'
  write ( *, * ) '  EXPONENTIAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  EXPONENTIAL_PDF evaluates the PDF.'

  x = 2.0

  a = 1.0
  b = 2.0

  call exponential_check ( a, b )

  call exponential_pdf ( x, a, b, pdf )

  call exponential_cdf ( x, a, b, cdf )

  call exponential_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value X           ', x2

  return
end
subroutine test057
!
!*******************************************************************************
!
!! TEST057 tests EXPONENTIAL_MEAN;
!! TEST057 tests EXPONENTIAL_SAMPLE;
!! TEST057 tests EXPONENTIAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST057'
  write ( *, * ) '  For the Exponential PDF:'
  write ( *, * ) '  EXPONENTIAL_MEAN computes the mean;'
  write ( *, * ) '  EXPONENTIAL_SAMPLE samples;'
  write ( *, * ) '  EXPONENTIAL_VARIANCE computes the variance.'

  a = 1.0
  b = 10.0

  call exponential_check ( a, b )

  call get_seed ( iseed )
  call exponential_mean ( a, b, mean )
  call exponential_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call exponential_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test059
!
!*******************************************************************************
!
!! TEST059 tests EXTREME_CDF.
!! TEST059 tests EXTREME_CDF_INV.
!! TEST059 tests EXTREME_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST059'
  write ( *, * ) '  For the Extreme CDF:'
  write ( *, * ) '  EXTREME_CDF evaluates the CDF;'
  write ( *, * ) '  EXTREME_CDF_INV inverts the CDF.'
  write ( *, * ) '  EXTREME_PDF evaluates the PDF;'

  x = 1.9

  a = 2.0
  b = 3.0

  call extreme_check ( a, b )

  call extreme_pdf ( x, a, b, pdf )

  call extreme_cdf ( x, a, b, cdf )

  call extreme_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test060
!
!*******************************************************************************
!
!! TEST060 tests EXTREME_MEAN;
!! TEST060 tests EXTREME_SAMPLE;
!! TEST060 tests EXTREME_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST060'
  write ( *, * ) '  For the Extreme PDF:'
  write ( *, * ) '  EXTREME_MEAN computes the mean;'
  write ( *, * ) '  EXTREME_SAMPLE samples;'
  write ( *, * ) '  EXTREME_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0

  call extreme_check ( a, b )

  call get_seed ( iseed )
  call extreme_mean ( a, b, mean )
  call extreme_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call extreme_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test062
!
!*******************************************************************************
!
!! TEST062 tests F_CENTRAL_CDF.
!! TEST062 tests F_CENTRAL_PDF.
!
  real cdf
  integer m
  integer n
  real pdf
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST062'
  write ( *, * ) '  For the central F PDF:'
  write ( *, * ) '  F_CENTRAL_CDF evaluates the CDF.'
  write ( *, * ) '  F_CENTRAL_PDF evaluates the PDF.'

  x = 648.0

  m = 1
  n = 1

  call f_central_pdf ( x, m, n, pdf )

  call f_central_cdf ( x, m, n, cdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =  ', x
  write ( *, * ) '  PDF parameter M = ', m
  write ( *, * ) '  PDF parameter N = ', n
  write ( *, * ) '  PDF value =       ', pdf
  write ( *, * ) '  CDF value =       ', cdf

  return
end
subroutine test063
!
!*******************************************************************************
!
!! TEST063 tests F_CENTRAL_MEAN;
!! TEST063 tests F_CENTRAL_SAMPLE;
!! TEST063 tests F_CENTRAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer iseed
  integer imax
  integer imin
  integer m
  real mean
  integer n
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST063'
  write ( *, * ) '  For the central F PDF:'
  write ( *, * ) '  F_CENTRAL_MEAN computes the mean;'
  write ( *, * ) '  F_CENTRAL_SAMPLE samples;'
  write ( *, * ) '  F_CENTRAL_VARIANCE computes the varianc.'

  m = 8
  n = 6

  call get_seed ( iseed )
  call f_central_mean ( m, n, mean )
  call f_central_variance ( m, n, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter M =             ', m
  write ( *, * ) '  PDF parameter N =             ', n
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call f_central_sample ( m, n, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test0645
!
!*******************************************************************************
!
!! TEST0645 tests FACTORIAL_LOG;
!! TEST0645 tests GAMMA_LOG_INT;
!
  real f
  real factorial_log
  real g
  real gamma_log_int
  integer i
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0645'
  write ( *, * ) '  FACTORIAL_LOG evaluates the log of the factorial function;'
  write ( *, * ) '  GAMMA_LOG_INT evaluates the log for integer argument.'

  write ( *, * ) ' '
  write ( *, * ) 'I, GAMMA_LOG_INT(I+1) FACTORIAL_LOG(I)'
  write ( *, * ) ' '

  do i = 1, 20
    g = gamma_log_int ( i+1 )
    f = factorial_log ( i )
    write ( *, '(i6,2g14.6)' ) i, g, f
  end do

  return
end
subroutine test065
!
!*******************************************************************************
!
!! TEST065 tests FACTORIAL_STIRLING;
!! TEST065 tests I_FACTORIAL;
!! TEST065 tests R_FACTORIAL.
!
  integer i_factorial
  real factorial_stirling
  integer i
  real r_factorial
  real value
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST065'
  write ( *, * ) '  FACTORIAL_STIRLING computes Stirling''s'
  write ( *, * ) '    approximate factorial function;'
  write ( *, * ) '  I_FACTORIAL evaluates the factorial function;'
  write ( *, * ) '  R_FACTORIAL evaluates the factorial function.'
  write ( *, * ) ' '
  write ( *, * ) '  N      Stirling     N!'
  write ( *, * ) ' '

  do i = 0, 10
    value = factorial_stirling ( i )
    write ( *, '(i6,g14.6,i20)' ) i, value, i_factorial ( i )
  end do

  write ( *, * ) ' '

  do i = 10, 20
    value = factorial_stirling ( i )
    write ( *, '(i6,2g14.6)' ) i, value, r_factorial ( i )
  end do

  return
end
subroutine test066
!
!*******************************************************************************
!
!! TEST066 tests FISK_CDF.
!! TEST066 tests FISK_CDF_INV.
!! TEST066 tests FISK_PDF.
!
  real a
  real b
  real c
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST066'
  write ( *, * ) '  For the Fisk PDF:'
  write ( *, * ) '  FISK_CDF evaluates the CDF;'
  write ( *, * ) '  FISK_CDF_INV inverts the CDF.'
  write ( *, * ) '  FISK_PDF evaluates the PDF;'

  x = 1.9

  a = 1.0
  b = 2.0
  c = 3.0

  call fisk_check ( a, b, c )

  call fisk_pdf ( x, a, b, c, pdf )

  call fisk_cdf ( x, a, b, c, cdf )

  call fisk_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test067
!
!*******************************************************************************
!
!! TEST067 tests FISK_MEAN;
!! TEST067 tests FISK_SAMPLE;
!! TEST067 tests FISK_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST067'
  write ( *, * ) '  For the Fisk PDF:'
  write ( *, * ) '  FISK_MEAN computes the mean;'
  write ( *, * ) '  FISK_SAMPLE samples;'
  write ( *, * ) '  FISK_VARIANCE computes the variance.'

  a = 1.0
  b = 2.0
  c = 3.0

  call fisk_check ( a, b, c )

  call get_seed ( iseed )
  call fisk_mean ( a, b, c, mean )
  call fisk_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call fisk_sample ( a, b, c, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test069
!
!*******************************************************************************
!
!! TEST069 tests FOLDED_NORMAL_CDF.
!! TEST069 tests FOLDED_NORMAL_CDF_INV.
!! TEST069 tests FOLDED_NORMAL_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST069'
  write ( *, * ) '  For the Folded Normal PDF:'
  write ( *, * ) '  FOLDED_NORMAL_CDF evaluates the CDF.'
  write ( *, * ) '  FOLDED_NORMAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  FOLDED_NORMAL_PDF evaluates the PDF.'

  x = 0.5

  a = 2.0
  b = 3.0

  call folded_normal_check ( a, b )

  call folded_normal_pdf ( x, a, b, pdf )

  call folded_normal_cdf ( x, a, b, cdf )

  call folded_normal_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value X =         ', x2

  return
end
subroutine test070
!
!*******************************************************************************
!
!! TEST070 tests FOLDED_NORMAL_MEAN;
!! TEST070 tests FOLDED_NORMAL_SAMPLE;
!! TEST070 tests FOLDED_NORMAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST070'
  write ( *, * ) '  For the Folded Normal PDF:'
  write ( *, * ) '  FOLDED_NORMAL_MEAN computes the mean;'
  write ( *, * ) '  FOLDED_NORMAL_SAMPLE samples;'
  write ( *, * ) '  FOLDED_NORMAL_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0

  call folded_normal_check ( a, b )

  call get_seed ( iseed )
  call folded_normal_mean ( a, b, mean )
  call folded_normal_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call folded_normal_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test072
!
!*******************************************************************************
!
!! TEST072 tests GAMMA;
!! TEST072 tests GAMMA_LOG;
!! TEST072 tests GAMMA_LOG_INT;
!! TEST072 tests R_FACTORIAL.
!
  real g1
  real g2
  real g3
  real g4
  real gamma
  real gamma_log
  real gamma_log_int
  integer i
  real r_factorial
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST072'
  write ( *, * ) '  GAMMA evaluates the Gamma function;'
  write ( *, * ) '  GAMMA_LOG evaluates the log of the Gamma function;'
  write ( *, * ) '  GAMMA_LOG_INT evaluates the log for integer argument;'
  write ( *, * ) '  R_FACTORIAL evaluates the factorial function.'

  write ( *, * ) ' '
  write ( *, * ) 'X, GAMMA(X), Exp(GAMMA_LOG(X)), Exp(GAMMA_LOG_INT(X)) ' // &
    'R_FACTORIAL(X+1)'
  write ( *, * ) ' '

  do i = 1, 10
    x = real ( i )
    g1 = gamma ( x )
    g2 = exp ( gamma_log ( x ) )
    g3 = exp ( gamma_log_int ( i ) )
    g4 = r_factorial ( i - 1 )
    write ( *, '(5g14.6)' ) x, g1, g2, g3, g4
  end do

  return
end
subroutine test073
!
!*******************************************************************************
!
!! TEST073 tests GAMMA_INC.
!
  integer, parameter :: test_num = 20
!
  real gamma_inc
  real p
  real p_test(test_num)
  integer test_i
  real x
  real x_test(test_num)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST073'
  write ( *, * ) '  GAMMA_INC evaluates the incomplete Gamma'
  write ( *, * ) '    function.'

  x_test( 1) =   2.0
  x_test( 2) =   6.0
  x_test( 3) =  10.0
  x_test( 4) =  20.0
  x_test( 5) =  20.0
  x_test( 6) =  50.0
  x_test( 7) =  50.0
  x_test( 8) =  50.0
  x_test( 9) =  50.0
  x_test(10) =   1.0
  x_test(11) =   5.0
  x_test(12) =   5.0
  x_test(13) =   7.0
  x_test(14) =  14.0
  x_test(15) =   2.0
  x_test(16) =  10.0
  x_test(17) =  20.0
  x_test(18) =  30.0
  x_test(19) =  93.0
  x_test(20) =   0.0

  p_test( 1) =   1.0
  p_test( 2) =   5.0
  p_test( 3) =   5.0
  p_test( 4) =   7.0
  p_test( 5) =  14.0
  p_test( 6) =   2.0
  p_test( 7) =  10.0
  p_test( 8) =  20.0
  p_test( 9) =  30.0
  p_test(10) =   2.0
  p_test(11) =   6.0
  p_test(12) =  10.0
  p_test(13) =  20.0
  p_test(14) =  20.0
  p_test(15) =  50.0
  p_test(16) =  50.0
  p_test(17) =  50.0
  p_test(18) =  50.0
  p_test(19) =  50.0
  p_test(20) =   1.0

  write ( *, * ) ' '
  write ( *, * ) '    X           P        GAMMA_INC(X,P)'
  write ( *, * ) ' '

  do test_i = 1, test_num

    x = x_test(test_i)
    p = p_test(test_i)
    write ( *, '(3g14.6)' ) x, p, gamma_inc ( x, p )
  
  end do

  return
end
subroutine test074
!
!*******************************************************************************
!
!! TEST074 tests GAMMA_CDF.
!! TEST074 tests GAMMA_PDF.
!
  real a
  real b
  real c
  real cdf
  integer i
  real pdf
  real x
!
  a = 1.0
  b = 1.5
  c = 3.0

  call gamma_check ( a, b, c )

  write ( *, * ) ' '
  write ( *, * ) 'TEST074'
  write ( *, * ) '  For the Gamma PDF:'
  write ( *, * ) '  GAMMA_CDF evaluates the CDF.'
  write ( *, * ) '  GAMMA_PDF evaluates the PDF.'
  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) ' '
  write ( *, * ) '  X  PDF   CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i ) / 5.0

    call gamma_cdf ( x, a, b, c, cdf )

    call gamma_pdf ( x, a, b, c, pdf )

    write ( *,     '(3g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test075
!
!*******************************************************************************
!
!! TEST075 tests GAMMA_MEAN;
!! TEST075 tests GAMMA_SAMPLE;
!! TEST075 tests GAMMA_VARIANCE.
!
  integer, parameter :: nsample = 1000
  integer, parameter :: test_num = 2
!
  real a
  real a_test(test_num)
  real b
  real b_test(test_num)
  real c
  real c_test(test_num)
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  integer test_i
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST075'
  write ( *, * ) '  For the Gamma PDF:'
  write ( *, * ) '  GAMMA_MEAN computes the mean;'
  write ( *, * ) '  GAMMA_SAMPLE samples;'
  write ( *, * ) '  GAMMA_VARIANCE computes the variance.'

  a_test(1) = 1.0
  a_test(2) = 2.0

  b_test(1) = 3.0
  b_test(2) = 0.5

  c_test(1) = 2.0
  c_test(2) = 0.5

  do test_i = 1, test_num

    a = a_test(test_i)
    b = b_test(test_i)
    c = c_test(test_i)

    call gamma_check ( a, b, c )
    call get_seed ( iseed )
    call gamma_mean ( a, b, c, mean )
    call gamma_variance ( a, b, c, variance )

    write ( *, * ) ' '
    write ( *, * ) '  Initial random number seed is ', iseed
    write ( *, * ) '  PDF parameter A =             ', a
    write ( *, * ) '  PDF parameter B =             ', b
    write ( *, * ) '  PDF parameter C =             ', c
    write ( *, * ) '  PDF mean =                    ', mean
    write ( *, * ) '  PDF variance =                ', variance
  
    do i = 1, nsample
      call gamma_sample ( a, b, c, iseed, x(i) )
    end do

    call rvec_mean ( nsample, x, mean )
    call rvec_variance ( nsample, x, variance )
    call rvec_max ( nsample, x, imax, xmax )
    call rvec_min ( nsample, x, imin, xmin )

    write ( *, * ) ' '
    write ( *, * ) '  Sample size =     ', nsample
    write ( *, * ) '  Sample mean =     ', mean
    write ( *, * ) '  Sample variance = ', variance
    write ( *, * ) '  Sample maximum =  ', xmax
    write ( *, * ) '  Sample minimum =  ', xmin

  end do

  return
end
subroutine test620
!
!*******************************************************************************
!
!! TEST620 tests GENLOGISTIC_CDF.
!! TEST620 tests GENLOGISTIC_CDF_INV.
!! TEST620 tests GENLOGISTIC_PDF.
!
  real a
  real b
  real c
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST620'
  write ( *, * ) '  For the Generalized Logistic PDF:'
  write ( *, * ) '  GENLOGISTIC_PDF evaluates the PDF.'
  write ( *, * ) '  GENLOGISTIC_CDF evaluates the CDF;'
  write ( *, * ) '  GENLOGISTIC_CDF_INV inverts the CDF.'

  x = 1.25

  a = 1.0
  b = 2.0
  c = 3.0

  call genlogistic_check ( a, b, c )

  call genlogistic_pdf ( x, a, b, c, pdf )

  call genlogistic_cdf ( x, a, b, c, cdf )

  call genlogistic_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value                     ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test621
!
!*******************************************************************************
!
!! TEST621 tests GENLOGISTIC_MEAN;
!! TEST621 tests GENLOGISTIC_SAMPLE;
!! TEST621 tests GENLOGISTIC_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST621'
  write ( *, * ) '  For the Generalized Logistic PDF:'
  write ( *, * ) '  GENLOGISTIC_MEAN computes the mean;'
  write ( *, * ) '  GENLOGISTIC_SAMPLE samples;'
  write ( *, * ) '  GENLOGISTIC_VARIANCE computes the variance.'

  a = 1.0
  b = 2.0
  c = 3.0

  call genlogistic_check ( a, b, c )

  call get_seed ( iseed )
  call genlogistic_mean ( a, b, c, mean )
  call genlogistic_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call genlogistic_sample ( a, b, c, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test077
!
!*******************************************************************************
!
!! TEST077 tests GEOMETRIC_CDF;
!! TEST077 tests GEOMETRIC_CDF_INV.
!! TEST077 tests GEOMETRIC_PDF;
!
  real a
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST077'
  write ( *, * ) '  For the Geometric PDF:'
  write ( *, * ) '  GEOMETRIC_CDF evaluates the CDF;'
  write ( *, * ) '  GEOMETRIC_CDF_INV inverts the CDF.'
  write ( *, * ) '  GEOMETRIC_PDF evaluates the PDF;'

  x = 5

  a = 0.25

  call geometric_check ( a )

  call geometric_pdf ( x, a, pdf )

  call geometric_cdf ( x, a, cdf )

  call geometric_cdf_inv ( cdf, a, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =       ', x
  write ( *, * ) '  PDF parameter A =      ', a
  write ( *, * ) '  PDF value =            ', pdf
  write ( *, * ) '  CDF value =            ', cdf
  write ( *, * ) '  CDF_INV value X =      ', x2

  return
end
subroutine test078
!
!*******************************************************************************
!
!! TEST078 tests GEOMETRIC_MEAN;
!! TEST078 tests GEOMETRIC_SAMPLE;
!! TEST078 tests GEOMETRIC_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST078'
  write ( *, * ) '  For the Geometric PDF:'
  write ( *, * ) '  GEOMETRIC_MEAN computes the mean;'
  write ( *, * ) '  GEOMETRIC_SAMPLE samples;'
  write ( *, * ) '  GEOMETRIC_VARIANCE computes the variance.'

  a = 0.25

  call geometric_check ( a )
  call get_seed ( iseed )
  call geometric_mean ( a, mean )
  call geometric_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call geometric_sample ( a, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test079
!
!*******************************************************************************
!
!! TEST079 tests GEOMETRIC_CDF.
!! TEST079 tests GEOMETRIC_PDF.
!
  real a
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST079'
  write ( *, * ) '  For the Geometric PDF:'
  write ( *, * ) '  GEOMETRIC_PDF evaluates the PDF.'
  write ( *, * ) '  GEOMETRIC_CDF evaluates the CDF.'

  a = 0.25

  call geometric_check ( a )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) ' '
  write ( *, * ) '      X      PDF(X)      CDF(X)'
  write ( *, * ) ' '

  do x = 0, 10
    call geometric_pdf ( x, a, pdf )
    call geometric_cdf ( x, a, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test0795
!
!*******************************************************************************
!
!! TEST0795 tests GOMPERTZ_CDF.
!! TEST0795 tests GOMPERTZ_CDF_INV.
!! TEST0795 tests GOMPERTZ_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0795'
  write ( *, * ) '  For the Gompertz PDF:'
  write ( *, * ) '  GOMPERTZ_CDF evaluates the CDF;'
  write ( *, * ) '  GOMPERTZ_CDF_INV inverts the CDF.'
  write ( *, * ) '  GOMPERTZ_PDF evaluates the PDF;'

  x = 0.6

  a = 2.0
  b = 3.0

  call gompertz_check ( a, b )

  call gompertz_pdf ( x, a, b, pdf )

  call gompertz_cdf ( x, a, b, cdf )

  call gompertz_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =        ', x
  write ( *, * ) '  PDF parameter A =       ', a
  write ( *, * ) '  PDF parameter B =       ', b
  write ( *, * ) '  PDF value =             ', pdf
  write ( *, * ) '  CDF value =             ', cdf
  write ( *, * ) '  CDF_INV value X =       ', x2

  return
end
subroutine test0796
!
!*******************************************************************************
!
!! TEST0796 tests GOMPERTZ_SAMPLE;
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0796'
  write ( *, * ) '  For the Gompertz PDF:'
  write ( *, * ) '  GOMPERTZ_SAMPLE samples;'

  a = 2.0
  b = 3.0

  call gompertz_check ( a, b )

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =       ', a
  write ( *, * ) '  PDF parameter B =       ', b
  
  do i = 1, nsample
    call gompertz_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test456
!
!*******************************************************************************
!
!! TEST456 tests GUMBEL_CDF;
!! TEST456 tests GUMBEL_CDF_INV.
!! TEST456 tests GUMBEL_PDF;
!
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST456'
  write ( *, * ) '  For the Gumbel PDF:'
  write ( *, * ) '  GUMBEL_CDF evaluates the CDF.'
  write ( *, * ) '  GUMBEL_CDF_INV inverts the CDF.'
  write ( *, * ) '  GUMBEL_PDF evaluates the PDF.'

  x = 0.5

  call gumbel_pdf ( x, pdf )

  call gumbel_cdf ( x, cdf )

  call gumbel_cdf_inv ( cdf, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value X           ', x2

  return
end
subroutine test457
!
!*******************************************************************************
!
!! TEST457 tests GUMBEL_MEAN;
!! TEST457 tests GUMBEL_SAMPLE;
!! TEST457 tests GUMBEL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST457'
  write ( *, * ) '  For the Gumbel PDF:'
  write ( *, * ) '  GUMBEL_MEAN computes the mean;'
  write ( *, * ) '  GUMBEL_SAMPLE samples;'
  write ( *, * ) '  GUMBEL_VARIANCE computes the variance.'

  call get_seed ( iseed )

  call gumbel_mean ( mean )

  call gumbel_variance ( variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF mean      =               ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call gumbel_sample ( iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test080
!
!*******************************************************************************
!
!! TEST080 tests HALF_NORMAL_CDF;
!! TEST080 tests HALF_NORMAL_CDF_INV.
!! TEST080 tests HALF_NORMAL_PDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST080'
  write ( *, * ) '  For the Half Normal PDF:'
  write ( *, * ) '  HALF_NORMAL_CDF evaluates the CDF.'
  write ( *, * ) '  HALF_NORMAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  HALF_NORMAL_PDF evaluates the PDF.'

  x = 0.5

  a = 0.0
  b = 2.0

  call half_normal_check ( a, b )

  call half_normal_pdf ( x, a, b, pdf )

  call half_normal_cdf ( x, a, b, cdf )

  call half_normal_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value X           ', x2

  return
end
subroutine test081
!
!*******************************************************************************
!
!! TEST081 tests HALF_NORMAL_MEAN;
!! TEST081 tests HALF_NORMAL_SAMPLE;
!! TEST081 tests HALF_NORMAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST081'
  write ( *, * ) '  For the Half Normal PDF:'
  write ( *, * ) '  HALF_NORMAL_MEAN computes the mean;'
  write ( *, * ) '  HALF_NORMAL_SAMPLE samples;'
  write ( *, * ) '  HALF_NORMAL_VARIANCE computes the variance.'

  a = 0.0
  b = 10.0

  call half_normal_check ( a, b )

  call get_seed ( iseed )
  call half_normal_mean ( a, b, mean )
  call half_normal_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call half_normal_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test083
!
!*******************************************************************************
!
!! TEST083 tests HYPERGEOMETRIC_CDF.
!! TEST083 tests HYPERGEOMETRIC_PDF.
!
  real cdf
  integer l
  integer m
  integer n
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST083'
  write ( *, * ) '  For the Hypergeometric PDF:'
  write ( *, * ) '  HYPERGEOMETRIC_CDF evaluates the CDF.'
  write ( *, * ) '  HYPERGEOMETRIC_PDF evaluates the PDF.'

  x = 7

  n = 100
  m = 70
  l = 1000

  call hypergeometric_check ( n, m, l )

  call hypergeometric_pdf ( x, n, m, l, pdf )

  call hypergeometric_cdf ( x, n, m, l, cdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =                ', x
  write ( *, * ) '  Total number of balls =         ', l
  write ( *, * ) '  Number of white balls =         ', m
  write ( *, * ) '  Number of balls taken =         ', n
  write ( *, * ) '  PDF value =                   = ', pdf
  write ( *, * ) '  CDF value =                   = ', cdf

  return
end
subroutine test085
!
!*******************************************************************************
!
!! TEST085 tests HYPERGEOMETRIC_MEAN;
!! TEST085 tests HYPERGEOMETRIC_SAMPLE;
!! TEST085 tests HYPERGEOMETRIC_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer iseed
  integer imax
  integer imin
  integer l
  integer m
  real mean
  integer n
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST085'
  write ( *, * ) '  For the Hypergeometric PDF:'
  write ( *, * ) '  HYPERGEOMETRIC_MEAN computes the mean;'
  write ( *, * ) '  HYPERGEOMETRIC_SAMPLE samples;'
  write ( *, * ) '  HYPERGEOMETRIC_VARIANCE computes the variance.'

  n = 100
  m = 70
  l = 1000

  call hypergeometric_check ( n, m, l )
  call get_seed ( iseed )
  call hypergeometric_mean ( n, m, l, mean )
  call hypergeometric_variance ( n, m, l, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter N =             ', n
  write ( *, * ) '  PDF parameter M =             ', m
  write ( *, * ) '  PDF parameter L =             ', l
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  write ( *, * ) ' '
  write ( *, * ) 'THIS CALL IS TAKING FOREVER!'
  return

  do i = 1, nsample
    call hypergeometric_sample ( n, m, l, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test086
!
!*******************************************************************************
!
!! TEST086 tests I_ROUNDUP.
!
  integer i
  integer i_roundup
  integer ival
  real rval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST086'
  write ( *, * ) '  I_ROUNDUP rounds reals up.'

  do i = -6, 6
    rval = real ( i ) / real ( 5.0 )
    ival = i_roundup ( rval )
    write ( *, '(g14.6,i6)' ) rval, ival
  end do

  return
end
subroutine test087
!
!*******************************************************************************
!
!! TEST087 tests INVERSE_GAUSSIAN_CDF.
!! TEST087 tests INVERSE_GAUSSIAN_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST087'
  write ( *, * ) '  For the Inverse Gaussian PDF:'
  write ( *, * ) '  INVERSE_GAUSSIAN_CDF evaluates the CDF.'
  write ( *, * ) '  INVERSE_GAUSSIAN_PDF evaluates the PDF.'

  x = 1.0

  a = 5.0
  b = 2.0

  call inverse_gaussian_check ( a, b )

  call inverse_gaussian_pdf ( x, a, b, pdf )

  call inverse_gaussian_cdf ( x, a, b, cdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X  = ', x
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) '  PDF parameter B = ', b
  write ( *, * ) '  PDF value =       ', pdf
  write ( *, * ) '  CDF value =       ', cdf

  return
end
subroutine test088
!
!*******************************************************************************
!
!! TEST088 tests INVERSE_GAUSSIAN_MEAN;
!! TEST088 tests INVERSE_GAUSSIAN_SAMPLE;
!! TEST088 tests INVERSE_GAUSSIAN_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST088'
  write ( *, * ) '  For the Inverse Gaussian PDF:'
  write ( *, * ) '  INVERSE_GAUSSIAN_MEAN computes the mean;'
  write ( *, * ) '  INVERSE_GAUSSIAN_SAMPLE samples;'
  write ( *, * ) '  INVERSE_GAUSSIAN_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0

  call inverse_gaussian_check ( a, b )

  call get_seed ( iseed )
  call inverse_gaussian_mean ( a, b, mean )
  call inverse_gaussian_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call inverse_gaussian_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test089
!
!*******************************************************************************
!
!! TEST089 tests LAPLACE_CDF.
!! TEST089 tests LAPLACE_CDF_INV.
!! TEST089 tests LAPLACE_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST089'
  write ( *, * ) '  For the Laplace PDF:'
  write ( *, * ) '  LAPLACE_CDF evaluates the CDF;'
  write ( *, * ) '  LAPLACE_CDF_INV inverts the CDF.'
  write ( *, * ) '  LAPLACE_PDF evaluates the PDF;'

  x = 3.0

  a = 1.0
  b = 2.0

  call laplace_check ( a, b )

  call laplace_pdf ( x, a, b, pdf )

  call laplace_cdf ( x, a, b, cdf )

  call laplace_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test090
!
!*******************************************************************************
!
!! TEST090 tests LAPLACE_MEAN;
!! TEST090 tests LAPLACE_SAMPLE;
!! TEST090 tests LAPLACE_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST090'
  write ( *, * ) '  For the Laplace PDF:'
  write ( *, * ) '  LAPLACE_MEAN computes the mean;'
  write ( *, * ) '  LAPLACE_SAMPLE samples;'
  write ( *, * ) '  LAPLACE_VARIANCE computes the variance.'

  a = 1.0
  b = 2.0

  call laplace_check ( a, b )

  call get_seed ( iseed )
  call laplace_mean ( a, b, mean )
  call laplace_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call laplace_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test092
!
!*******************************************************************************
!
!! TEST092 tests LOGISTIC_CDF.
!! TEST092 tests LOGISTIC_CDF_INV.
!! TEST092 tests LOGISTIC_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST092'
  write ( *, * ) '  For the Logistic PDF:'
  write ( *, * ) '  LOGISTIC_CDF evaluates the CDF;'
  write ( *, * ) '  LOGISTIC_CDF_INV inverts the CDF.'
  write ( *, * ) '  LOGISTIC_PDF evaluates the PDF;'

  x = 3.0

  a = 1.0
  b = 2.0

  call logistic_check ( a, b )

  call logistic_cdf ( x, a, b, cdf )

  call logistic_pdf ( x, a, b, pdf )

  call logistic_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test093
!
!*******************************************************************************
!
!! TEST093 tests LOGISTIC_MEAN;
!! TEST093 tests LOGISTIC_SAMPLE;
!! TEST093 tests LOGISTIC_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST093'
  write ( *, * ) '  For the Logistic PDF:'
  write ( *, * ) '  LOGISTIC_MEAN computes the mean;'
  write ( *, * ) '  LOGISTIC_SAMPLE samples;'
  write ( *, * ) '  LOGISTIC_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0

  call logistic_check ( a, b )

  call get_seed ( iseed )
  call logistic_mean ( a, b, mean )
  call logistic_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call logistic_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test095
!
!*******************************************************************************
!
!! TEST095 tests LOGNORMAL_CDF.
!! TEST095 tests LOGNORMAL_CDF_INV.
!! TEST095 tests LOGNORMAL_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST095'
  write ( *, * ) '  For the Lognormal PDF:'
  write ( *, * ) '  LOGNORMAL_CDF evaluates the CDF;'
  write ( *, * ) '  LOGNORMAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  LOGNORMAL_PDF evaluates the PDF;'

  x = exp ( 9.0 )

  a = 10.0
  b = 2.25

  call lognormal_check ( a, b )

  call lognormal_cdf ( x, a, b, cdf )

  call lognormal_pdf ( x, a, b, pdf )

  call lognormal_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test096
!
!*******************************************************************************
!
!! TEST096 tests LOGNORMAL_MEAN;
!! TEST096 tests LOGNORMAL_SAMPLE;
!! TEST096 tests LOGNORMAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST096'
  write ( *, * ) '  For the Lognormal PDF:'
  write ( *, * ) '  LOGNORMAL_MEAN computes the mean;'
  write ( *, * ) '  LOGNORMAL_SAMPLE samples;'
  write ( *, * ) '  LOGNORMAL_VARIANCE computes the variance.'

  a = 1.0
  b = 2.0

  call lognormal_check ( a, b )

  call get_seed ( iseed )
  call lognormal_mean ( a, b, mean )
  call lognormal_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call lognormal_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test098
!
!*******************************************************************************
!
!! TEST098 tests LOGSERIES_CDF.
!! TEST098 tests LOGSERIES_CDF_INV.
!! TEST098 tests LOGSERIES_PDF.
!
  real a
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST098'
  write ( *, * ) '  For the Logseries PDF,'
  write ( *, * ) '  LOGSERIES_CDF evaluates the CDF;'
  write ( *, * ) '  LOGSERIES_CDF_INV inverts the CDF.'
  write ( *, * ) '  LOGSERIES_PDF evaluates the PDF;'

  x = 3

  a = 0.25

  call logseries_check ( a )

  call logseries_pdf ( x, a, pdf )

  call logseries_cdf ( x, a, cdf )

  call logseries_cdf_inv ( cdf, a, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument =     ', x
  write ( *, * ) '  PDF parameter A =  ', a
  write ( *, * ) '  PDF value =        ', pdf
  write ( *, * ) '  CDF value =        ', cdf
  write ( *, * ) '  CDF_INV value   =  ', x2

  return
end
subroutine test099
!
!*******************************************************************************
!
!! TEST099 tests LOGSERIES_CDF.
!! TEST099 tests LOGSERIES_PDF.
!
  real a
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST099'
  write ( *, * ) '  For the Logseries PDF:'
  write ( *, * ) '  LOGSERIES_CDF evaluates the CDF;'
  write ( *, * ) '  LOGSERIES_PDF evaluates the PDF.'

  x = 2

  a = 0.25

  call logseries_check ( a )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X  = ', x
  write ( *, * ) ' '
  write ( *, * ) '  X        PDF(X)       CDF(X)'
  write ( *, * ) ' '

  do x = 1, 10
    call logseries_pdf ( x, a, pdf )
    call logseries_cdf ( x, a, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test100
!
!*******************************************************************************
!
!! TEST100 tests LOGSERIES_MEAN.
!! TEST100 tests LOGSERIES_SAMPLE.
!! TEST100 tests LOGSERIES_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST100'
  write ( *, * ) '  For the Logseries PDF:'
  write ( *, * ) '  LOGSERIES_MEAN computes the mean;'
  write ( *, * ) '  LOGSERIES_VARIANCE computes the variance;'
  write ( *, * ) '  LOGSERIES_SAMPLE samples.'

  a = 0.25

  call logseries_check ( a )
  call logseries_mean ( a, mean )
  call logseries_variance ( a, variance )

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call logseries_sample ( a, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test101
!
!*******************************************************************************
!
!! TEST101 tests LORENTZ_CDF.
!! TEST101 tests LORENTZ_CDF_INV.
!! TEST101 tests LORENTZ_PDF.
!
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST101'
  write ( *, * ) '  For the Lorentz PDF:'
  write ( *, * ) '  LORENTZ_CDF evaluates the CDF;'
  write ( *, * ) '  LORENTZ_CDF_INV inverts the CDF.'
  write ( *, * ) '  LORENTZ_PDF evaluates the PDF;'

  x = 0.75

  call lorentz_pdf ( x, pdf )

  call lorentz_cdf ( x, cdf )

  call lorentz_cdf_inv ( cdf, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF value =               ', pdf
  write ( *, * ) '  CDF value =               ', cdf
  write ( *, * ) '  CDF_INV value X =         ', x2

  return
end
subroutine test102
!
!*******************************************************************************
!
!! TEST102 tests LORENTZ_MEAN;
!! TEST102 tests LORENTZ_SAMPLE;
!! TEST102 tests LORENTZ_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST102'
  write ( *, * ) '  For the Lorentz PDF:'
  write ( *, * ) '  LORENTZ_MEAN computes the mean;'
  write ( *, * ) '  LORENTZ_VARIANCE computes the variance;'
  write ( *, * ) '  LORENTZ_SAMPLE samples.'

  call get_seed ( iseed )
  call lorentz_mean ( mean )
  call lorentz_variance ( variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call lorentz_sample ( iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test104
!
!*******************************************************************************
!
!! TEST104 tests MAXWELL_CDF.
!! TEST104 tests MAXWELL_CDF_INV.
!! TEST104 tests MAXWELL_PDF.
!
  real a
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST104'
  write ( *, * ) '  For the Maxwell CDF:'
  write ( *, * ) '  MAXWELL_CDF evaluates the CDF.'
  write ( *, * ) '  MAXWELL_CDF_INV inverts the CDF.'
  write ( *, * ) '  MAXWELL_PDF evaluates the PDF.'

  x = 0.5

  a = 2.0

  call maxwell_check ( a )

  call maxwell_pdf ( x, a, pdf )

  call maxwell_cdf ( x, a, cdf )

  call maxwell_cdf_inv ( cdf, a, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value             ', x2

  return
end
subroutine test105
!
!*******************************************************************************
!
!! TEST105 tests MAXWELL_MEAN;
!! TEST105 tests MAXWELL_SAMPLE;
!! TEST105 tests MAXWELL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST105'
  write ( *, * ) '  For the Maxwell PDF:'
  write ( *, * ) '  MAXWELL_MEAN computes the mean;'
  write ( *, * ) '  MAXWELL_VARIANCE computes the variance;'
  write ( *, * ) '  MAXWELL_SAMPLE samples.'

  a = 2.0

  call get_seed ( iseed )
  call maxwell_mean ( a, mean )
  call maxwell_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF mean =                    ', variance

  do i = 1, nsample
    call maxwell_sample ( a, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test106
!
!*******************************************************************************
!
!! TEST106 tests MULTINOMIAL_COEF1;
!! TEST106 tests MULTINOMIAL_COEF2.
!
  integer, parameter :: maxfactor = 5
!
  integer factor(maxfactor)
  integer i
  integer j
  integer n
  integer ncomb1
  integer ncomb2
  integer nfactor
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST106'
  write ( *, * ) '  MULTINOMIAL_COEF1 computes multinomial'
  write ( *, * ) '    coefficients using the Gamma function;'
  write ( *, * ) '  MULTINOMIAL_COEF2 computes multinomial'
  write ( *, * ) '    coefficients directly.'

  write ( *, * ) ' '
  write ( *, * ) '  Line 10 of the BINOMIAL table:'
  write ( *, * ) ' '

  n = 10
  nfactor = 2

  do i = 0, n

    factor(1) = i
    factor(2) = n - i

    call multinomial_coef1 ( nfactor, factor, ncomb1 )

    call multinomial_coef2 ( nfactor, factor, ncomb2 )

    write ( *, '(i4,i4,3x,i5,i5)' ) factor(1), factor(2), ncomb1, ncomb2

  end do

  write ( *, * ) ' '
  write ( *, * ) '  Level 5 of the TRINOMIAL coefficients:'

  n = 5
  nfactor = 3

  do i = 0, n

    factor(1) = i

    write ( *, * ) ' '

    do j = 0, n - factor(1)

      factor(2) = j
      factor(3) = n - factor(1) - factor(2)

      call multinomial_coef1 ( nfactor, factor, ncomb1 )

      call multinomial_coef2 ( nfactor, factor, ncomb2 )

      write ( *, '(i4,i4,i4,3x,i5,i5)' ) factor(1), factor(2), factor(3), &
        ncomb1, ncomb2

    end do

  end do

  return
end
subroutine test107
!
!*******************************************************************************
!
!! TEST107 tests MULTINOMIAL_MEAN;
!! TEST107 tests MULTINOMIAL_SAMPLE;
!! TEST107 tests MULTINOMIAL_VARIANCE;
!! TEST107 tests IROW_MEAN;
!! TEST107 tests IROW_VARIANCE.
!
  integer, parameter :: b = 3
  integer, parameter :: nsample = 1000
!
  integer a
  real c(b)
  integer i
  integer iseed
  integer imax(b)
  integer imin(b)
  real mean(b)
  real variance(b)
  integer x(b,nsample)
  integer xmax(b)
  integer xmin(b)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST107'
  write ( *, * ) '  For the Multinomial PDF:'
  write ( *, * ) '  MULTINOMIAL_MEAN computes the mean;'
  write ( *, * ) '  MULTINOMIAL_SAMPLE samples;'
  write ( *, * ) '  MULTINOMIAL_VARIANCE computes the variance;'

  a = 5

  c(1) = 0.125
  c(2) = 0.500
  c(3) = 0.375

  call multinomial_check ( a, b, c )

  call get_seed ( iseed )
  call multinomial_mean ( a, b, c, mean )
  call multinomial_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C = '
  do i = 1, b
    write ( *, '(g14.6)' ) c(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  PDF means and variances:'
  write ( *, * ) ' '
  do i = 1, b
    write ( *, '(2g14.6)' ) mean(i), variance(i)
  end do
  
  do i = 1, nsample
    call multinomial_sample ( a, b, c, iseed, x(1,i) )
  end do

  call irow_max ( b, b, nsample, x, imax, xmax )
  call irow_min ( b, b, nsample, x, imin, xmin )
  call irow_mean ( b, b, nsample, x, mean )
  call irow_variance ( b, b, nsample, x, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size = ', nsample
  write ( *, * ) '  Component Min, Max, Mean, Variance:'
  do i = 1, b
    write ( *, '(3i6,2g14.6)' ) i, xmin(i), xmax(i), mean(i), variance(i)
  end do

  return
end
subroutine test108
!
!*******************************************************************************
!
!! TEST108 tests MULTINOMIAL_PDF.
!
  integer, parameter :: b = 3
!
  integer a
  real c(b)
  integer i
  real pdf
  integer x(b)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST108'
  write ( *, * ) '  For the Multinomial PDF:'
  write ( *, * ) '  MULTINOMIAL_PDF evaluates the PDF.'

  x(1) = 0
  x(2) = 2
  x(3) = 3

  a = 5

  c(1) = 0.10
  c(2) = 0.50
  c(3) = 0.40

  call multinomial_check ( a, b, c )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X:'
  write ( *, * ) ' '
  do i = 1, b
    write ( *, '(i6)' ) x(i)
  end do
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C = '
  do i = 1, b
    write ( *, '(g14.6)' ) c(i)
  end do

  call multinomial_pdf ( x, a, b, c, pdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF value =     ', pdf

  return
end
subroutine test520
!
!*******************************************************************************
!
!! TEST520 tests NAKAGAMI_CDF.
!! TEST520 tests NAKAGAMI_PDF.
!
  real a
  real b
  real c
  real cdf
  real pdf
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST520'
  write ( *, * ) '  For the Nakagami PDF:'
  write ( *, * ) '  NAKAGAMI_CDF evaluates the CDF;'
  write ( *, * ) '  NAKAGAMI_PDF evaluates the PDF;'

  x = 1.25

  a = 1.0
  b = 2.0
  c = 3.0

  call nakagami_check ( a, b, c )

  call nakagami_pdf ( x, a, b, c, pdf )

  call nakagami_cdf ( x, a, b, c, cdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf

  return
end
subroutine test521
!
!*******************************************************************************
!
!! TEST521 tests NAKAGAMI_MEAN;
!! TEST521 tests NAKAGAMI_VARIANCE.
!
  real a
  real b
  real c
  integer iseed
  real mean
  real variance
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST521'
  write ( *, * ) '  For the Nakagami PDF:'
  write ( *, * ) '  NAKAGAMI_MEAN computes the mean;'
  write ( *, * ) '  NAKAGAMI_VARIANCE computes the variance.'

  a = 1.0
  b = 2.0
  c = 3.0

  call nakagami_check ( a, b, c )

  call get_seed ( iseed )
  call nakagami_mean ( a, b, c, mean )
  call nakagami_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  return
end
subroutine test109
!
!*******************************************************************************
!
!! TEST109 tests NORMAL_01_CDF;
!! TEST109 tests NORMAL_01_CDF_INV.
!! TEST109 tests NORMAL_01_PDF;
!
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST109'
  write ( *, * ) '  For the Normal 01 PDF:'
  write ( *, * ) '  NORMAL_01_CDF evaluates the CDF;'
  write ( *, * ) '  NORMAL_01_CDF_INV inverts the CDF.'
  write ( *, * ) '  NORMAL_01_PDF evaluates the PDF;'

  x = 1.0

  call normal_01_pdf ( x, pdf )

  call normal_01_cdf ( x, cdf )

  call normal_01_cdf_inv ( cdf, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =       ', x
  write ( *, * ) '  PDF value =            ', pdf
  write ( *, * ) '  CDF value =            ', cdf
  write ( *, * ) '  CDF_INV value X =      ', x2

  return
end
subroutine test110
!
!*******************************************************************************
!
!! TEST110 tests NORMAL_01_MEAN;
!! TEST110 tests NORMAL_01_SAMPLE;
!! TEST110 tests NORMAL_01_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST110'
  write ( *, * ) '  For the Normal 01 PDF:'
  write ( *, * ) '  NORMAL_01_MEAN computes the mean;'
  write ( *, * ) '  NORMAL_01_SAMPLE samples the PDF;'
  write ( *, * ) '  NORMAL_01_VARIANCE returns the variance.'

  call get_seed ( iseed )
  call normal_01_mean ( mean )
  call normal_01_variance ( variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call normal_01_sample ( iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test112
!
!*******************************************************************************
!
!! TEST112 tests NORMAL_CDF;
!! TEST112 tests NORMAL_CDF_INV.
!! TEST112 tests NORMAL_PDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST112'
  write ( *, * ) '  For the Normal PDF:'
  write ( *, * ) '  NORMAL_CDF evaluates the CDF;'
  write ( *, * ) '  NORMAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  NORMAL_PDF evaluates the PDF;'

  x = 90.0

  a = 100.0
  b = 15.0

  call normal_check ( a, b )

  call normal_pdf ( x, a, b, pdf )

  call normal_cdf ( x, a, b, cdf )

  call normal_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =       ', x
  write ( *, * ) '  PDF parameter A =      ', a
  write ( *, * ) '  PDF parameter B =      ', b
  write ( *, * ) '  PDF value =            ', pdf
  write ( *, * ) '  CDF value =            ', cdf
  write ( *, * ) '  CDF_INV value X =      ', x2

  return
end
subroutine test113
!
!*******************************************************************************
!
!! TEST113 tests NORMAL_MEAN;
!! TEST113 tests NORMAL_SAMPLE;
!! TEST113 tests NORMAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST113'
  write ( *, * ) '  For the Normal PDF:'
  write ( *, * ) '  NORMAL_MEAN computes the mean;'
  write ( *, * ) '  NORMAL_SAMPLE samples;'
  write ( *, * ) '  NORMAL_VARIANCE returns the variance.'

  a = 100.0
  b = 15.0

  call normal_check ( a, b )

  call get_seed ( iseed )
  call normal_mean ( a, b, mean )
  call normal_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call normal_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test115
!
!*******************************************************************************
!
!! TEST115 tests PARETO_CDF.
!! TEST115 tests PARETO_CDF_INV.
!! TEST115 tests PARETO_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST115'
  write ( *, * ) '  For the Pareto PDF:'
  write ( *, * ) '  PARETO_CDF evaluates the CDF;'
  write ( *, * ) '  PARETO_CDF_INV inverts the CDF.'
  write ( *, * ) '  PARETO_PDF evaluates the PDF;'

  x = 4.0

  a = 2.0
  b = 3.0

  call pareto_check ( a, b )

  call pareto_pdf ( x, a, b, pdf )

  call pareto_cdf ( x, a, b, cdf )

  call pareto_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test116
!
!*******************************************************************************
!
!! TEST116 tests PARETO_MEAN;
!! TEST116 tests PARETO_SAMPLE;
!! TEST116 tests PARETO_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST116'
  write ( *, * ) '  For the Pareto PDF:'
  write ( *, * ) '  PARETO_MEAN computes the mean;'
  write ( *, * ) '  PARETO_SAMPLE samples;'
  write ( *, * ) '  PARETO_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0

  call pareto_check ( a, b )

  call get_seed ( iseed )
  call pareto_mean ( a, b, mean )
  call pareto_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call pareto_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test118
!
!*******************************************************************************
!
!! TEST118 tests PASCAL_MEAN;
!! TEST118 tests PASCAL_SAMPLE;
!! TEST118 tests PASCAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST118'
  write ( *, * ) '  For the Pascal PDF:'
  write ( *, * ) '  PASCAL_MEAN computes the mean;'
  write ( *, * ) '  PASCAL_SAMPLE samples;'
  write ( *, * ) '  PASCAL_VARIANCE computes the variance.'

  a = 2
  b = 0.75

  call pascal_check ( a, b )

  call get_seed ( iseed )
  call pascal_mean ( a, b, mean )
  call pascal_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call pascal_sample ( a, b, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test119
!
!*******************************************************************************
!
!! TEST119 tests PASCAL_CDF.
!! TEST119 tests PASCAL_CDF_INV.
!! TEST119 tests PASCAL_PDF.
!
  integer a
  real b
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST119'
  write ( *, * ) '  For the Pascal PDF:'
  write ( *, * ) '  PASCAL_CDF evaluates the CDF.'
  write ( *, * ) '  PASCAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  PASCAL_PDF evaluates the PDF.'

  a = 2
  b = 0.25

  call pascal_check ( a, b )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) '  PDF parameter B = ', b
  write ( *, * ) ' '
  write ( *, * ) '      X      PDF(X)       CDF(X)    CDF_INV(CDF)'
  write ( *, * ) ' '

  do x = 0, 10
    call pascal_pdf ( x, a, b, pdf )
    call pascal_cdf ( x, a, b, cdf )
    call pascal_cdf_inv ( cdf, a, b, x2 )
    write ( *, '(i6,2g14.6,i6)' ) x, pdf, cdf, x2
  end do

  return
end
subroutine test120
!
!*******************************************************************************
!
!! TEST120 tests PEARSON_05_PDF.
!
  real a
  real b
  real c
  real pdf
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST120'
  write ( *, * ) '  For the Pearson 05 PDF:'
  write ( *, * ) '  PEARSON_05_PDF evaluates the PDF.'

  x = 5.0

  a = 1.0
  b = 2.0
  c = 3.0

  call pearson_05_check ( a, b, c )
 
  call pearson_05_pdf ( x, a, b, c, pdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =  ', x
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) '  PDF parameter B = ', b
  write ( *, * ) '  PDF parameter C = ', c
  write ( *, * ) '  PDF value =       ', pdf

  return
end
subroutine test121
!
!*******************************************************************************
!
!! TEST121 tests PLANCK_PDF.
!
  real pdf
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST121'
  write ( *, * ) '  For the Planck PDF:'
  write ( *, * ) '  PLANCK_PDF evaluates the PDF.'

  x = 0.5

  call planck_pdf ( x, pdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =  ', x
  write ( *, * ) '  PDF value =       ', pdf

  return
end
subroutine test122
!
!*******************************************************************************
!
!! TEST122 tests PLANCK_SAMPLE.
!
  integer, parameter :: nsample = 1000
!
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST122'
  write ( *, * ) '  For the Planck PDF:'
  write ( *, * ) '  PLANCK_SAMPLE samples.'

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  
  do i = 1, nsample
    call planck_sample ( iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test123
!
!*******************************************************************************
!
!! TEST123 tests POISSON_CDF.
!! TEST123 tests POISSON_CDF_INV.
!
  real a
  real cdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST123'
  write ( *, * ) '  For the Poisson PDF:'
  write ( *, * ) '  POISSON_CDF evaluates the CDF,'
  write ( *, * ) '  POISSON_CDF_INV inverts the CDF.'

  x = 7

  a = 10.0

  call poisson_check ( a )

  call poisson_cdf ( x, a, cdf )

  call poisson_cdf_inv ( cdf, a, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =   ', x
  write ( *, * ) '  PDF parameter A =  ', a
  write ( *, * ) '  CDF value =        ', cdf
  write ( *, * ) '  CDF_INV value X =  ', x2

  return
end
subroutine test124
!
!*******************************************************************************
!
!! TEST124 tests POISSON_MEAN;
!! TEST124 tests POISSON_SAMPLE;
!! TEST124 tests POISSON_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST124'
  write ( *, * ) '  For the Poisson PDF:'
  write ( *, * ) '  POISSON_MEAN computes the mean;'
  write ( *, * ) '  POISSON_SAMPLE samples;'
  write ( *, * ) '  POISSON_VARIANCE computes the variance.'

  a = 10.0

  call poisson_check ( a )
  call get_seed ( iseed )
  call poisson_mean ( a, mean )
  call poisson_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call poisson_sample ( a, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test125
!
!*******************************************************************************
!
!! TEST125 tests POISSON_CDF.
!! TEST125 tests POISSON_PDF.
!
  real a
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST125'
  write ( *, * ) '  For the Poisson PDF:'
  write ( *, * ) '  POISSON_PDF evaluates the PDF.'
  write ( *, * ) '  POISSON_CDF evaluates the CDF.'

  a = 10.0

  call poisson_check ( a )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) ' '
  write ( *, * ) '      X      PDF(X)      CDF(X)'
  write ( *, * ) ' '

  do x = 0, 10
    call poisson_pdf ( x, a, pdf )
    call poisson_cdf ( x, a, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test126
!
!*******************************************************************************
!
!! TEST126 tests POLY_VAL.
!
  integer, parameter :: n = 3
!
  real coef(0:n)
  integer i
  real val
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST126'
  write ( *, * ) '  POLY_VAL evaluates a polynomial.'

  coef(0) = 7.0
  coef(1) = 5.0
  coef(2) = 2.0
  coef(3) = 1.0

  write ( *, * ) ' '
  write ( *, * ) '  Polynomial degree is N = ', n
  write ( *, * ) '  Polynomial coefficients are:'
  do i = 0, n
    write ( *, '(i6,g14.6)' ) i, coef(i)
  end do

  x = 2.0

  call poly_val ( n, coef, x, val )

  write ( *, * ) ' '
  write ( *, * ) '  Polynomial argument X = ', x
  write ( *, * ) '  Polynomial value VAL =  ', val
  write ( *, * ) ' '
  write ( *, * ) '  (Expected value is 33.)'

  return
end
subroutine test209
!
!*******************************************************************************
!
!! TEST209 tests POWER_CDF;
!! TEST209 tests POWER_CDF_INV.
!! TEST209 tests POWER_PDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST209'
  write ( *, * ) '  For the Power PDF:'
  write ( *, * ) '  POWER_CDF evaluates the CDF;'
  write ( *, * ) '  POWER_CDF_INV inverts the CDF.'
  write ( *, * ) '  POWER_PDF evaluates the PDF;'

  x = 0.6

  a = 2.0
  b = 3.0

  call power_check ( a, b )

  call power_pdf ( x, a, b, pdf )

  call power_cdf ( x, a, b, cdf )

  call power_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =        ', x
  write ( *, * ) '  PDF parameter A =       ', a
  write ( *, * ) '  PDF parameter B =       ', b
  write ( *, * ) '  PDF value =             ', pdf
  write ( *, * ) '  CDF value =             ', cdf
  write ( *, * ) '  CDF_INV value X =       ', x2

  return
end
subroutine test211
!
!*******************************************************************************
!
!! TEST211 tests POWER_MEAN;
!! TEST211 tests POWER_SAMPLE;
!! TEST211 tests POWER_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST211'
  write ( *, * ) '  For the Power PDF:'
  write ( *, * ) '  POWER_MEAN computes the mean;'
  write ( *, * ) '  POWER_SAMPLE samples;'
  write ( *, * ) '  POWER_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0

  call power_check ( a, b )

  call get_seed ( iseed )
  call power_mean ( a, b, mean )
  call power_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call power_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test127
!
!*******************************************************************************
!
!! TEST127 tests RAYLEIGH_CDF.
!! TEST127 tests RAYLEIGH_CDF_INV.
!! TEST127 tests RAYLEIGH_PDF.
!
  real a
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST127'
  write ( *, * ) '  For the Rayleigh PDF:'
  write ( *, * ) '  RAYLEIGH_CDF evaluates the CDF;'
  write ( *, * ) '  RAYLEIGH_CDF_INV inverts the CDF.'
  write ( *, * ) '  RAYLEIGH_PDF evaluates the PDF;'

  x = 1.0

  a = 2.0

  call rayleigh_check ( a )

  call rayleigh_pdf ( x, a, pdf )

  call rayleigh_cdf ( x, a, cdf )

  call rayleigh_cdf_inv ( cdf, a, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test128
!
!*******************************************************************************
!
!! TEST128 tests RAYLEIGH_MEAN;
!! TEST128 tests RAYLEIGH_SAMPLE;
!! TEST128 tests RAYLEIGH_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST128'
  write ( *, * ) '  For the Rayleigh PDF:'
  write ( *, * ) '  RAYLEIGH_MEAN computes the mean;'
  write ( *, * ) '  RAYLEIGH_SAMPLE samples;'
  write ( *, * ) '  RAYLEIGH_VARIANCE computes the variance.'

  a = 2.0

  call rayleigh_check ( a )
  call get_seed ( iseed )
  call rayleigh_mean ( a, mean )
  call rayleigh_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call rayleigh_sample ( a, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test380
!
!*******************************************************************************
!
!! TEST380 tests RECIPROCAL_CDF;
!! TEST380 tests RECIPROCAL_CDF_INV.
!! TEST380 tests RECIPROCAL_CDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST380'
  write ( *, * ) '  For the Reciprocal PDF:'
  write ( *, * ) '  RECIPROCAL_CDF evaluates the CDF.'
  write ( *, * ) '  RECIPROCAL_CDF_INV inverts the CDF.'
  write ( *, * ) '  RECIPROCAL_PDF evaluates the PDF.'

  x = 1.5

  a = 1.0
  b = 3.0

  call reciprocal_check ( a, b )

  call reciprocal_pdf ( x, a, b, pdf )

  call reciprocal_cdf ( x, a, b, cdf )

  call reciprocal_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value X           ', x2

  return
end
subroutine test381
!
!*******************************************************************************
!
!! TEST381 tests RECIPROCAL_MEAN;
!! TEST381 tests RECIPROCAL_SAMPLE;
!! TEST381 tests RECIPROCAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST381'
  write ( *, * ) '  For the Reciprocal PDF:'
  write ( *, * ) '  RECIPROCAL_MEAN computes the mean;'
  write ( *, * ) '  RECIPROCAL_SAMPLE samples;'
  write ( *, * ) '  RECIPROCAL_VARIANCE computes the variance.'

  a = 1.0
  b = 3.0

  call reciprocal_check ( a, b )

  call get_seed ( iseed )
  call reciprocal_mean ( a, b, mean )
  call reciprocal_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call reciprocal_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test304
!
!*******************************************************************************
!
!! TEST304 tests SECH_CDF.
!! TEST304 tests SECH_CDF_INV.
!! TEST304 tests SECH_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST304'
  write ( *, * ) '  For the Hyperbolic Secant PDF:'
  write ( *, * ) '  SECH_CDF evaluates the CDF.'
  write ( *, * ) '  SECH_CDF_INV inverts the CDF.'
  write ( *, * ) '  SECH_PDF evaluates the PDF.'

  x = 4.0

  a = 3.0
  b = 2.0

  call sech_check ( a, b )

  call sech_pdf ( x, a, b, pdf )

  call sech_cdf ( x, a, b, cdf )

  call sech_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDF_INV value =           ', x2

  return
end
subroutine test305
!
!*******************************************************************************
!
!! TEST305 tests SECH_MEAN;
!! TEST305 tests SECH_SAMPLE;
!! TEST305 tests SECH_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST305'
  write ( *, * ) '  For the Hyperbolic Secant PDF:'
  write ( *, * ) '  SECH_MEAN computes the mean;'
  write ( *, * ) '  SECH_SAMPLE samples;'
  write ( *, * ) '  SECH_VARIANCE computes the variance.'

  a = 3.0
  b = 2.0

  call sech_check ( a, b )

  call get_seed ( iseed )
  call sech_mean ( a, b, mean )
  call sech_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call sech_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test204
!
!*******************************************************************************
!
!! TEST204 tests SEMICIRCULAR_CDF.
!! TEST204 tests SEMICIRCULAR_CDF_INV.
!! TEST204 tests SEMICIRCULAR_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST204'
  write ( *, * ) '  For the Semicircular PDF:'
  write ( *, * ) '  SEMICIRCULAR_CDF evaluates the CDF.'
  write ( *, * ) '  SEMICIRCULAR_CDF_INV inverts the CDF.'
  write ( *, * ) '  SEMICIRCULAR_PDF evaluates the PDF.'

  x = 4.0

  a = 3.0
  b = 2.0

  call semicircular_check ( a, b )

  call semicircular_pdf ( x, a, b, pdf )

  call semicircular_cdf ( x, a, b, cdf )

  call semicircular_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =          ', x
  write ( *, * ) '  PDF parameter A =         ', a
  write ( *, * ) '  PDF parameter B =         ', b
  write ( *, * ) '  PDF value                 ', pdf
  write ( *, * ) '  CDF value                 ', cdf
  write ( *, * ) '  CDV_INV value X =         ', x2

  return
end
subroutine test205
!
!*******************************************************************************
!
!! TEST205 tests SEMICIRCULAR_MEAN;
!! TEST205 tests SEMICIRCULAR_SAMPLE;
!! TEST205 tests SEMICIRCULAR_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST205'
  write ( *, * ) '  For the Semicircular PDF:'
  write ( *, * ) '  SEMICIRCULAR_MEAN computes the mean;'
  write ( *, * ) '  SEMICIRCULAR_SAMPLE samples;'
  write ( *, * ) '  SEMICIRCULAR_VARIANCE computes the variance.'

  a = 3.0
  b = 2.0

  call semicircular_check ( a, b )

  call get_seed ( iseed )
  call semicircular_mean ( a, b, mean )
  call semicircular_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call semicircular_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test130
!
!*******************************************************************************
!
!! TEST130 tests STUDENT_CENTRAL_CDF.
!! TEST130 tests STUDENT_CENTRAL_PDF.
!
  real a
  real b
  real c
  real cdf
  real pdf
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST130'
  write ( *, * ) '  For the central Student PDF:'
  write ( *, * ) '  STUDENT_CENTRAL_CDF evaluates the CDF.'
  write ( *, * ) '  STUDENT_CENTRAL_PDF evaluates the PDF.'

  x = 2.447

  a = 0.5
  b = 2.0
  c = 6.0

  call student_central_check ( a, b, c )

  call student_central_pdf ( x, a, b, c, pdf )

  call student_central_cdf ( x, a, b, c, cdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =    ', x
  write ( *, * ) '  PDF parameter A =   ', a
  write ( *, * ) '  PDF parameter B =   ', b
  write ( *, * ) '  PDF parameter C =   ', c
  write ( *, * ) '  PDF value =         ', pdf
  write ( *, * ) '  CDF value =         ', cdf

  return
end
subroutine test131
!
!*******************************************************************************
!
!! TEST131 tests STUDENT_CENTRAL_MEAN;
!! TEST131 tests STUDENT_CENTRAL_SAMPLE;
!! TEST131 tests STUDENT_CENTRAL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST131'
  write ( *, * ) '  For the central Student PDF:'
  write ( *, * ) '  STUDENT_CENTRAL_MEAN computes the mean;'
  write ( *, * ) '  STUDENT_CENTRAL_SAMPLE samples;'
  write ( *, * ) '  STUDENT_CENTRAL_VARIANCE computes the variance.'

  a = 0.5
  b = 2.0
  c = 6.0

  call student_central_check ( a, b, c )
  call get_seed ( iseed )
  call student_central_mean ( a, b, c, mean )
  call student_central_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance
  
  do i = 1, nsample
    call student_central_sample ( a, b, c, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test133
!
!*******************************************************************************
!
!! TEST133 tests STUDENT_NONCENTRAL_CDF.
!
  real b
  real cdf
  integer idf
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST133'
  write ( *, * ) '  For the Noncentral Student PDF:'
  write ( *, * ) '  STUDENT_NONCENTRAL_CDF evaluates the CDF;'

  x = 0.50

  idf = 10
  b = 1.0

  call student_noncentral_cdf ( x, idf, b, cdf )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter IDF =           ', idf
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  CDF value =                   ', cdf

  return
end
subroutine test134
!
!*******************************************************************************
!
!! TEST134 tests TRIANGULAR_CDF;
!! TEST134 tests TRIANGULAR_CDF_INV.
!! TEST134 tests TRIANGULAR_PDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST134'
  write ( *, * ) '  For the Triangular PDF:'
  write ( *, * ) '  TRIANGULAR_CDF evaluates the CDF;'
  write ( *, * ) '  TRIANGULAR_CDF_INV inverts the CDF.'
  write ( *, * ) '  TRIANGULAR_PDF evaluates the PDF;'

  x = 4.0

  a = 1.0
  b = 10.0

  call triangular_check ( a, b )

  call triangular_pdf ( x, a, b, pdf )

  call triangular_cdf ( x, a, b, cdf )

  call triangular_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =       ', x
  write ( *, * ) '  PDF parameter A =      ', a
  write ( *, * ) '  PDF parameter B =      ', b
  write ( *, * ) '  PDF value =            ', pdf
  write ( *, * ) '  CDF value =            ', cdf
  write ( *, * ) '  CDF_INV value X =      ', x2

  return
end
subroutine test135
!
!*******************************************************************************
!
!! TEST135 tests TRIANGULAR_MEAN;
!! TEST135 tests TRIANGULAR_SAMPLE;
!! TEST135 tests TRIANGULAR_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST135'
  write ( *, * ) '  For the Triangular PDF:'
  write ( *, * ) '  TRIANGULAR_MEAN computes mean;'
  write ( *, * ) '  TRIANGULAR_SAMPLE samples;'
  write ( *, * ) '  TRIANGULAR_VARIANCE computes variance.'

  a = 1.0
  b = 10.0

  call triangular_check ( a, b )

  call get_seed ( iseed )
  call triangular_mean ( a, b, mean )
  call triangular_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call triangular_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test137
!
!*******************************************************************************
!
!! TEST137 tests UNIFORM_01_ORDER_SAMPLE;
!
  integer, parameter :: n = 10
!
  integer i
  integer iseed
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST137'
  write ( *, * ) '  For the Uniform 01 Order PDF:'
  write ( *, * ) '  UNIFORM_ORDER_SAMPLE samples.'

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed

  call uniform_01_order_sample ( n, iseed, x )

  write ( *, * ) ' '
  write ( *, * ) 'Ordered sample:'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, x(i)
  end do

  return
end
subroutine test138
!
!*******************************************************************************
!
!! TEST138 tests UNIFORM_NSPHERE_SAMPLE;
!
  integer, parameter :: n = 3
!
  integer i
  integer iseed
  integer j
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST138'
  write ( *, * ) '  For the Uniform PDF on the N-Sphere:'
  write ( *, * ) '  UNIFORM_NSPHERE_SAMPLE samples.'

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  Dimension N of sphere =       ', n
  write ( *, * ) ' '
  write ( *, * ) 'Points on the sphere:'
  write ( *, * ) ' '

  do i = 1, 10
    call uniform_nsphere_sample ( n, iseed, x )
    write ( *, '(i6,3g14.6)' ) i, x(1:n)
  end do

  return
end
subroutine test139
!
!*******************************************************************************
!
!! TEST139 tests UNIFORM_CDF;
!! TEST139 tests UNIFORM_CDF_INV.
!! TEST139 tests UNIFORM_PDF;
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST139'
  write ( *, * ) '  For the Uniform PDF:'
  write ( *, * ) '  UNIFORM_CDF evaluates the CDF;'
  write ( *, * ) '  UNIFORM_CDF_INV inverts the CDF.'
  write ( *, * ) '  UNIFORM_PDF evaluates the PDF;'

  x = 4.0

  a = 1.0
  b = 10.0

  call uniform_check ( a, b )

  call uniform_pdf ( x, a, b, pdf )

  call uniform_cdf ( x, a, b, cdf )

  call uniform_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =       ', x
  write ( *, * ) '  PDF parameter A =      ', a
  write ( *, * ) '  PDF parameter B =      ', b
  write ( *, * ) '  PDF value =            ', pdf
  write ( *, * ) '  CDF value =            ', cdf
  write ( *, * ) '  CDF_INV value X =      ', x2

  return
end
subroutine test140
!
!*******************************************************************************
!
!! TEST140 tests UNIFORM_MEAN;
!! TEST140 tests UNIFORM_SAMPLE;
!! TEST140 tests UNIFORM_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST140'
  write ( *, * ) '  For the Uniform PDF:'
  write ( *, * ) '  UNIFORM_MEAN computes mean;'
  write ( *, * ) '  UNIFORM_SAMPLE samples;'
  write ( *, * ) '  UNIFORM_VARIANCE computes variance.'

  a = 1.0
  b = 10.0

  call uniform_check ( a, b )
  call get_seed ( iseed )
  call uniform_mean ( a, b, mean )
  call uniform_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call uniform_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test142
!
!*******************************************************************************
!
!! TEST142 tests UNIFORM_DISCRETE_CDF;
!! TEST142 tests UNIFORM_DISCRETE_CDF_INV.
!! TEST142 tests UNIFORM_DISCRETE_PDF;
!
  integer a
  integer b
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST142'
  write ( *, * ) '  For the Uniform Discrete PDF:'
  write ( *, * ) '  UNIFORM_DISCRETE_CDF evaluates the CDF;'
  write ( *, * ) '  UNIFORM_DISCRETE_CDF_INV inverts the CDF.'
  write ( *, * ) '  UNIFORM_DISCRETE_PDF evaluates the PDF;'

  x = 4

  a = 1
  b = 6

  call uniform_discrete_check ( a, b )

  call uniform_discrete_pdf ( x, a, b, pdf )

  call uniform_discrete_cdf ( x, a, b, cdf )

  call uniform_discrete_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X2 =            ', x2

  return
end
subroutine test143
!
!*******************************************************************************
!
!! TEST143 tests UNIFORM_DISCRETE_MEAN;
!! TEST143 tests UNIFORM_DISCRETE_SAMPLE;
!! TEST143 tests UNIFORM_DISCRETE_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  integer a
  integer b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST143'
  write ( *, * ) '  For the Uniform discrete PDF:'
  write ( *, * ) '  UNIFORM_DISCRETE_MEAN computes the mean;'
  write ( *, * ) '  UNIFORM_DISCRETE_SAMPLE samples;'
  write ( *, * ) '  UNIFORM_DISCRETE_VARIANCE computes the variance.'

  a = 1
  b = 6

  call uniform_discrete_check ( a, b )
  call get_seed ( iseed )
  call uniform_discrete_mean ( a, b, mean )
  call uniform_discrete_variance ( a, b, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call uniform_discrete_sample ( a, b, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test144
!
!*******************************************************************************
!
!! TEST144 tests UNIFORM_DISCRETE_CDF.
!! TEST144 tests UNIFORM_DISCRETE_PDF.
!
  integer a
  integer b
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST144'
  write ( *, * ) '  For the Uniform discrete PDF.'
  write ( *, * ) '  UNIFORM_DISCRETE_PDF evaluates the PDF.'
  write ( *, * ) '  UNIFORM_DISCRETE_CDF evaluates the CDF.'

  a = 1
  b = 6

  call uniform_discrete_check ( a, b )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) ' '
  write ( *, * ) '      X      PDF(X)      CDF(X)'
  write ( *, * ) ' '

  do x = 0, 6
    call uniform_discrete_pdf ( x, a, b, pdf )
    call uniform_discrete_cdf ( x, a, b, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test145
!
!*******************************************************************************
!
!! TEST145 tests VON_MISES_CDF.
!! TEST145 tests VON_MISES_CDF_INV.
!! TEST145 tests VON_MISES_PDF.
!
  real a
  real b
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST145'
  write ( *, * ) '  For the Von Mises PDF:'
  write ( *, * ) '  VON_MISES_CDF evaluates the CDF.'
  write ( *, * ) '  VON_MISES_CDF_INV inverts the CDF.'
  write ( *, * ) '  VON_MISES_PDF evaluates the PDF.'

  x = 0.5

  a = 1.0
  b = 2.0

  call von_mises_check ( a, b )

  call von_mises_pdf ( x, a, b, pdf )

  call von_mises_cdf ( x, a, b, cdf )

  call von_mises_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =       ', x
  write ( *, * ) '  PDF parameter A =      ', a
  write ( *, * ) '  PDF parameter B =      ', b
  write ( *, * ) '  PDF value =            ', pdf
  write ( *, * ) '  CDF value =            ', cdf
  write ( *, * ) '  CDF_INV value X =      ', x2

  return
end
subroutine test146
!
!*******************************************************************************
!
!! TEST146 tests VON_MISES_MEAN;
!! TEST146 tests VON_MISES_SAMPLE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST146'
  write ( *, * ) '  For the Von Mises PDF:'
  write ( *, * ) '  VON_MISES_MEAN computes the mean;'
  write ( *, * ) '  VON_MISES_SAMPLE samples.'

  a = 1.0
  b = 2.0

  call von_mises_check ( a, b )
  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b

  do i = 1, nsample
    call von_mises_sample ( a, b, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test148
!
!*******************************************************************************
!
!! TEST148 tests WEIBULL_CDF.
!! TEST148 tests WEIBULL_CDF_INV.
!! TEST148 tests WEIBULL_PDF.
!
  real a
  real b
  real c
  real cdf
  real pdf
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST148'
  write ( *, * ) '  For the Weibull PDF:'
  write ( *, * ) '  WEIBULL_CDF evaluates the CDF;'
  write ( *, * ) '  WEIBULL_CDF_INV inverts the CDF.'
  write ( *, * ) '  WEIBULL_PDF evaluates the PDF;'

  x = 3.0

  a = 2.0
  b = 3.0
  c = 4.0

  call weibull_check ( a, b, c )

  call weibull_pdf ( x, a, b, c, pdf )

  call weibull_cdf ( x, a, b, c, cdf )

  call weibull_cdf_inv ( cdf, a, b, c, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test149
!
!*******************************************************************************
!
!! TEST149 tests WEIBULL_MEAN;
!! TEST149 tests WEIBULL_SAMPLE;
!! TEST149 tests WEIBULL_VARIANCE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  real c
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  real x(nsample)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST149'
  write ( *, * ) '  For the Weibull PDF:'
  write ( *, * ) '  WEIBULL_MEAN computes the mean;'
  write ( *, * ) '  WEIBULL_SAMPLE samples;'
  write ( *, * ) '  WEIBULL_VARIANCE computes the variance.'

  a = 2.0
  b = 3.0
  c = 4.0

  call weibull_check ( a, b, c )

  call get_seed ( iseed )
  call weibull_mean ( a, b, c, mean )
  call weibull_variance ( a, b, c, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF parameter C =             ', c
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call weibull_sample ( a, b, c, iseed, x(i) )
  end do

  call rvec_mean ( nsample, x, mean )
  call rvec_variance ( nsample, x, variance )
  call rvec_max ( nsample, x, imax, xmax )
  call rvec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test151
!
!*******************************************************************************
!
!! TEST151 tests WEIBULL_DISCRETE_CDF.
!! TEST151 tests WEIBULL_DISCRETE_CDF_INV.
!! TEST151 tests WEIBULL_DISCRETE_PDF.
!
  real a
  real b
  real cdf
  real pdf
  integer x
  integer x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST151'
  write ( *, * ) '  For the Weibull Discrete PDF,'
  write ( *, * ) '  WEIBULL_DISCRETE_CDF evaluates the CDF;'
  write ( *, * ) '  WEIBULL_DISCRETE_CDF_INV inverts the CDF.'
  write ( *, * ) '  WEIBULL_DISCRETE_PDF evaluates the PDF;'

  x = 2

  a = 0.50
  b = 1.5

  call weibull_discrete_check ( a, b )

  call weibull_discrete_pdf ( x, a, b, pdf )

  call weibull_discrete_cdf ( x, a, b, cdf )

  call weibull_discrete_cdf_inv ( cdf, a, b, x2 )

  write ( *, * ) ' '
  write ( *, * ) '  PDF argument X =              ', x
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF parameter B =             ', b
  write ( *, * ) '  PDF value =                   ', pdf
  write ( *, * ) '  CDF value =                   ', cdf
  write ( *, * ) '  CDF_INV value X =             ', x2

  return
end
subroutine test152
!
!*******************************************************************************
!
!! TEST152 tests WEIBULL_DISCRETE_CDF.
!! TEST152 tests WEIBULL_DISCRETE_PDF.
!
  real a
  real b
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST152'
  write ( *, * ) '  For the Weibull Discrete PDF:'
  write ( *, * ) '  WEIBULL_DISCRETE_PDF evaluates the PDF;'
  write ( *, * ) '  WEIBULL_DISCRETE_CDF evaluates the CDF.'

  a = 0.50
  b = 1.5

  call weibull_discrete_check ( a, b )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A =     ', a
  write ( *, * ) '  PDF parameter B =     ', b
  write ( *, * ) ' '
  write ( *, * ) '      X      PDF(X)      CDF(X)'
  write ( *, * ) ' '

  do x = 0, 10
    call weibull_discrete_pdf ( x, a, b, pdf )
    call weibull_discrete_cdf ( x, a, b, cdf )
    write ( *, '(i6,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test153
!
!*******************************************************************************
!
!! TEST153 tests WEIBULL_DISCRETE_SAMPLE.
!
  integer, parameter :: nsample = 1000
!
  real a
  real b
  integer i
  integer imax
  integer imin
  integer iseed
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST153'
  write ( *, * ) '  For the discrete Weibull PDF:'
  write ( *, * ) '  WEIBULL_DISCRETE_SAMPLE samples.'

  a = 0.5
  b = 1.5

  call weibull_discrete_check ( a, b )

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =     ', a
  write ( *, * ) '  PDF parameter B =     ', b
  
  do i = 1, nsample
    call weibull_discrete_sample ( a, b, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
subroutine test154
!
!*******************************************************************************
!
!! TEST154 tests ZIPF_CDF.
!! TEST154 tests ZIPF_PDF.
!
  real a
  real cdf
  real pdf
  integer x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST154'
  write ( *, * ) '  For the Zipf PDF:'
  write ( *, * ) '  ZIPF_PDF evaluates the PDF.'
  write ( *, * ) '  ZIPF_CDF evaluates the CDF.'

  a = 2.0

  call zipf_check ( a )

  write ( *, * ) ' '
  write ( *, * ) '  PDF parameter A = ', a
  write ( *, * ) ' '
  write ( *, * ) '  X    PDF(X)       CDF(X)'
  write ( *, * ) ' '

  do x = 1, 20

    call zipf_pdf ( x, a, pdf )
    call zipf_cdf ( x, a, cdf )
    write ( *, '(i6,2x,2g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test155
!
!*******************************************************************************
!
!! TEST155 tests ZIPF_SAMPLE.
!
  integer, parameter :: nsample = 1000
!
  real a
  integer i
  integer iseed
  integer imax
  integer imin
  real mean
  real variance
  integer x(nsample)
  integer xmax
  integer xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST155'
  write ( *, * ) '  For the Zipf PDF:'
  write ( *, * ) '  ZIPF_SAMPLE samples.'

  a = 4.0

  call zipf_check ( a )
  call get_seed ( iseed )
  call zipf_mean ( a, mean )
  call zipf_variance ( a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Initial random number seed is ', iseed
  write ( *, * ) '  PDF parameter A =             ', a
  write ( *, * ) '  PDF mean =                    ', mean
  write ( *, * ) '  PDF variance =                ', variance

  do i = 1, nsample
    call zipf_sample ( a, iseed, x(i) )
  end do

  call ivec_mean ( nsample, x, mean )
  call ivec_variance ( nsample, x, variance )
  call ivec_max ( nsample, x, imax, xmax )
  call ivec_min ( nsample, x, imin, xmin )

  write ( *, * ) ' '
  write ( *, * ) '  Sample size =     ', nsample
  write ( *, * ) '  Sample mean =     ', mean
  write ( *, * ) '  Sample variance = ', variance
  write ( *, * ) '  Sample maximum =  ', xmax
  write ( *, * ) '  Sample minimum =  ', xmin

  return
end
