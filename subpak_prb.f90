program subpak_prb
!
!*******************************************************************************
!
!! SUBPAK_PRB tests routines from the SUBPAK library.
!
  character ( len = 8 ) date
  character ( len = 10 ) time
!
  call date_and_time ( date, time )

  write ( *, * ) ' '
  write ( *, * ) 'SUBPAK_PRB'
  write ( *, * ) '  A set of test programs for SUBPAK.'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
 
  call test001
  call test002
  call test003
  call test004
  call test005
  call test006
  call test007
  call test008
  call test0083
  call test0085

  call test011
  call test012
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
  call test023
  call test024
  call test025
  call test026
  call test027
  call test028
  call test029
  call test030

  call test031
  call test032
  call test033
  call test888
  call test034
  call test0345
  call test035
  call test036
  call test0365
  call test0366
  call test037
  call test038
  call test0385
  call test039
  call test040

  call test041
  call test042
  call test043
  call test0435
  call test044
  call test0445
  call test045
  call test0451
  call test0452
  call test0453
  call test046
  call test047
  call test0475
  call test048
  call test049
  call test050

  call test051
  call test052
  call test053
  call test054
  call test055
  call test056
  call test057
  call test058
  call test059
  call test0595
  call test0596
  call test05965
  call test0597
  call test060

  call test061
  call test0615
  call test062
  call test063
  call test0635
  call test064
  call test065
  call test066
  call test067
  call test0675
  call test0676
  call test0677
  call test068
  call test1007
  call test0685
  call test069
  call test070

  call test071
  call test072
  call test073
  call test1735
  call test0775
  call test078
  call test0745
  call test075
  call test0756
  call test009
  call test0754
  call test0755
  call test076
  call test0761
  call test082
  call test083
  call test084
  call test0734
  call test0735
  call test079
  call test074
  call test0799
  call test07995
  call test080
  call test077

  call test081
  call test0813
  call test0814
  call test0815
  call test085
  call test086
  call test087
  call test088
  call test089
  call test090

  call test091
  call test092
  call test093
  call test094
  call test095
  call test096
  call test097
  call test098
  call test0985
  call test099
  call test100
  call test1005

  call test101
  call test102
  call test103
  call test1035
  call test104
  call test105
  call test106
  call test107
  call test1071
  call test1072
  call test1073
  call test108
  call test109
  call test110

  call test111
  call test112
  call test113
  call test114
  call test115
  call test1154
  call test1155
  call test116
  call test1165
  call test126
  call test117
  call test1174
  call test1175
  call test118
  call test119
  call test120

  call test121
  call test122
  call test123
  call test124
  call test1245
  call test127
  call test128

  write ( *, * ) ' '
  write ( *, * ) 'SUBPAK_PRB'
  write ( *, * ) '  Normal end of SUBPAK tests.'

  stop
end
subroutine test001
!
!*******************************************************************************
!
!! TEST001 tests ACOSH2.
!
  real a
  real acosh2
  integer i
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST001'
  write ( *, * ) '  ACOSH2 computes the inverse hyperbolic cosine'
  write ( *, * ) '  of a given value.'
  write ( *, * ) ' '
  write ( *, * ) '  X     ACOSH2(X)     COSH(ACOSH2(X))'
  write ( *, * ) ' '

  do i = 0, 10
    x = 1.0 + real ( i ) / 5.0
    a = acosh2 ( x )
    x2 = cosh ( a )
    write ( *, '(3g14.6)' ) x, a, x2
  end do

  return
end
subroutine test002
!
!*******************************************************************************
!
!! TEST002 tests AGUD.
!! TEST002 tests GUD.
!
  real agud
  real gamma
  real gud
  integer i
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST002'
  write ( *, * ) '  AGUD computes the inverse Gudermannian;'
  write ( *, * ) '  GUD computes the Gudermannian.'
  write ( *, * ) ' '
  write ( *, * ) '  X     GUD(X)     AGUD(GUD(X))'
  write ( *, * ) ' '

  do i = 0, 10
    x = 1.0 + real ( i ) / 5.0
    gamma = gud ( x )
    x2 = agud ( gamma )
    write ( *, '(3g14.6)' ) x, gamma, x2
  end do

  return
end
subroutine test003
!
!*******************************************************************************
!
!! TEST003 tests ASINH2.
!
  real a
  real asinh2
  integer i
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST003'
  write ( *, * ) '  ASINH2 computes the inverse hyperbolic sine'
  write ( *, * ) '  of a given value.'
  write ( *, * ) ' '
  write ( *, * ) '  X     ASINH2(X)     SINH(ASINH2(X))'
  write ( *, * ) ' '

  do i = 0, 10
    x = 1.0 + real ( i ) / 5.0
    a = asinh2 ( x )
    x2 = sinh ( a )
    write ( *, '(3g14.6)' ) x, a, x2
  end do

  return
end
subroutine test004
!
!*******************************************************************************
!
!! TEST004 tests ATAN4.
!
  integer, parameter :: ntest = 8
!
  real atan4
  integer i
  real x
  real xtest(ntest)
  real y
  real ytest(ntest)
!
  xtest(1) = 1.0E+00
  ytest(1) = 0.0E+00

  xtest(2) = 1.0E+00
  ytest(2) = 1.0E+00

  xtest(3) = 0.0E+00
  ytest(3) = 1.0E+00

  xtest(4) = -1.0E+00
  ytest(4) = 1.0E+00

  xtest(5) = -1.0E+00
  ytest(5) = 0.0E+00

  xtest(6) = - 1.0E+00
  ytest(6) = - 1.0E+00

  xtest(7) =   0.0E+00
  ytest(7) = - 1.0E+00

  xtest(8) =   1.0E+00
  ytest(8) = - 1.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST004'
  write ( *, * ) '  ATAN4 computes the arc-tangent given Y and X;'
  write ( *, * ) '  ATAN2 is the system version of this routine.'
  write ( *, * ) ' '
  write ( *, * ) '  X     Y     ATAN2(Y,X)   ATAN4(Y,X)'
  write ( *, * ) ' '

  do i = 1, ntest
    x = xtest(i)
    y = ytest(i)
    write ( *, '(4g14.6)' ) x, y, atan2 ( y, x ), atan4 ( y, x )
  end do

  return
end
subroutine test005
!
!*******************************************************************************
!
!! TEST005 tests ATANH2.
!
  real a
  real atanh2
  integer i
  real x
  real x2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST005'
  write ( *, * ) '  ATANH2 computes the inverse hyperbolic tangent'
  write ( *, * ) '    of a given value.'
  write ( *, * ) ' '
  write ( *, * ) '  X     ATANH2(X)     TANH(ATANH2(X))'
  write ( *, * ) ' '

  do i = -2, 9
    x = real ( i ) / 10.0E+00
    a = atanh2 ( x )
    x2 = tanh ( a )
    write ( *, '(3g14.6)' ) x, a, x2
  end do

  return
end
subroutine test006
!
!*******************************************************************************
!
!! TEST006 tests AXIS_LIMITS.
!
  integer ndivs
  integer nticks
  real pxdiv
  real pxmax
  real pxmin
  real xmax
  real xmin
!
  xmin = 67.3E+00
  xmax = 114.7E+00
  ndivs = 6

  call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

  write ( *, * ) ' '
  write ( *, * ) 'TEST006'
  write ( *, * ) '  AXIS_LIMITS adjusts plot limits'
  write ( *, * ) '    to "nicer" values.'
  write ( *, * ) ' '
  write ( *, * ) '  Input XMIN, XMAX = ', xmin, xmax
  write ( *, * ) '  Input NDIVS = ', ndivs
  write ( *, * ) ' '
  write ( *, * ) '  Output PXMIN, PXMAX = ', pxmin, pxmax
  write ( *, * ) '  Output PXDIV = ', pxdiv
  write ( *, * ) '  Output NTICKS = ', nticks

  xmin = - 26.0E+00
  xmax = + 26.0E+00
  ndivs = 10

  call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

  write ( *, * ) ' '
  write ( *, * ) '  Input XMIN, XMAX = ', xmin, xmax
  write ( *, * ) '  Input NDIVS = ', ndivs
  write ( *, * ) ' '
  write ( *, * ) '  Output PXMIN, PXMAX = ', pxmin, pxmax
  write ( *, * ) '  Output PXDIV = ', pxdiv
  write ( *, * ) '  Output NTICKS = ', nticks

  return
end
subroutine test007
!
!*******************************************************************************
!
!! TEST007 tests AXIS_LIMITS.
!
  integer, parameter :: ntest = 5
!
  integer i
  integer ndivs
  integer nticks
  real pxdiv
  real pxmax
  real pxmin
  real test_max(ntest)
  real test_min(ntest)
  real xmax
  real xmin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST007'
  write ( *, * ) '  AXIS_LIMITS computes "nice" limits for a graph'
  write ( *, * ) '    that must include a given range.'

  ndivs = 5

  test_min(1) = 1.0E+00
  test_max(1) = 9.0E+00

  test_min(2) = 1.003E+00
  test_max(2) = 4.125E+00

  test_min(3) = 101.25E+00
  test_max(3) = 193.75E+00

  test_min(4) = 2000.125E+00
  test_max(4) = 2000.250E+00

  test_min(5) = -7.0E+00
  test_max(5) = 12.0E+00

  write ( *, * ) ' '
  write ( *, * ) '  All tests use NDIVS = ', ndivs
  write ( *, * ) ' '
  write ( *, * ) '  XMIN  XMAX  PXMIN  PXMAX  PXDIV NTICKS'
  write ( *, * ) ' '

  do i = 1, ntest

    xmin = test_min(i)
    xmax = test_max(i)

    call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, nticks )

    write ( *, '(5g12.4,i6)' ) xmin, xmax, pxmin, pxmax, pxdiv, nticks

  end do

  return
end
subroutine test008
!
!*******************************************************************************
!
!! TEST008 tests BAR_CHECK;
!! TEST008 tests BAR_CODE;
!! TEST008 tests BAR_DIGIT_CODE.
!
  character ( len = 113 ) bar
  integer digit(12)
  character ( len = 7 ) codel
  character ( len = 7 ) coder
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST008'
  write ( *, * ) '  BAR_CHECK checks digits for a barcode;'
  write ( *, * ) '  BAR_CODE computes the barcode for a string of'
  write ( *, * ) '    11 digits;'
  write ( *, * ) '  BAR_DIGIT_CODE returns the left and right codes'
  write ( *, * ) '    for each digit.'

  do i = 1, 11
    digit(i) = mod ( i-1, 10 )
  end do
 
  call bar_check ( digit )
 
  write ( *, * ) ' '
  write ( *, * ) '  Check digit is ', digit(12)
 
  write ( *, * ) ' '
  write ( *, * ) '  Digit code:'
  write ( *, * ) ' '
  do i = 0, 9
    call bar_digit_code ( i, codel, coder )
    write ( *, '(i2,2x,a7,2x,a7)' ) i, codel, coder
  end do
 
  call bar_code ( digit, bar )
 
  write ( *, * ) ' '
  write ( *, * ) '  Bar code:'
  write ( *, * ) ' '
  write ( *, '(a)' ) bar(1:9)
  write ( *, '(a)' ) bar(10:12)
  write ( *, '(a)' ) bar(13:19)
  write ( *, '(a)' ) bar(20:26)
  write ( *, '(a)' ) bar(27:33)
  write ( *, '(a)' ) bar(34:40)
  write ( *, '(a)' ) bar(41:47)
  write ( *, '(a)' ) bar(48:54)
  write ( *, '(a)' ) bar(55:59)
  write ( *, '(a)' ) bar(60:66)
  write ( *, '(a)' ) bar(67:73)
  write ( *, '(a)' ) bar(74:80)
  write ( *, '(a)' ) bar(81:87)
  write ( *, '(a)' ) bar(88:94)
  write ( *, '(a)' ) bar(95:101)
  write ( *, '(a)' ) bar(102:104)
  write ( *, '(a)' ) bar(105:113)
 
  return
end
subroutine test0083
!
!*******************************************************************************
!
!! TEST0083 tests BMI_ENGLISH.
!
  real bmi
  real bmi_english
  real h
  real h_ft
  real h_in
  integer i
  real w
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0083'
  write ( *, * ) '  BMI_ENGLISH computes the Body Mass Index'
  write ( *, * ) '  given body measurements in English Units.'
  write ( *, * ) ' '
  write ( *, * ) '  Weight(LB)  Height (FT/IN)  BMI'
  write ( *, * ) ' '

  do i = 1, 10

    call r_random ( 100.0E+00, 250.0E+00, w )
    call r_random ( 4.0E+00, 6.75E+00, h )

    h_ft = int ( h )
    h_in = nint ( 12.0E+00 * ( h - h_ft ) )
 
    bmi = bmi_english ( w, h_ft, h_in )
    write ( *, '(4f10.2)' ) w, h_ft, h_in, bmi

  end do

  return
end
subroutine test0085
!
!*******************************************************************************
!
!! TEST0085 tests CHVEC_PERMUTE.
!
  integer, parameter :: n = 10
!
  character chvec(n)
  integer i
  integer iseed
  integer p(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0085'
  write ( *, * ) '  CHVEC_PERMUTE applies a permutation to a character vector.'

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Using random number seed ISEED = ', iseed

  call perm_random ( n, p )

  call ivec_print ( n, p, '  The random permutation:' )

  do i = 1, n
    chvec(i) = char ( ichar ( 'A' ) + i - 1 )
  end do

  call chvec_print ( n, chvec, '  CHVEC before permutation:' )

  call chvec_permute ( n, chvec, p )

  call chvec_print ( n, chvec, '  CHVEC after permutation:' )
 
  return
end
subroutine test011
!
!*******************************************************************************
!
!! TEST011 tests FAC_DIV;
!! TEST011 tests FAC_MUL;
!! TEST011 tests FAC_LCM;
!! TEST011 tests FAC_GCD;
!! TEST011 tests I_TO_FAC.
!
  integer, parameter :: nprime = 5
!
  integer i1
  integer i2
  integer ibot
  integer itop
  integer npower1(nprime)
  integer npower2(nprime)
  integer npower3(nprime)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST011'
  write ( *, * ) '  For products of prime factors:'
  write ( *, * ) '  FAC_DIV computes a quotient;'
  write ( *, * ) '  FAC_MUL multiplies;'
  write ( *, * ) '  FAC_LCM computes the LCM;'
  write ( *, * ) '  FAC_GCD computes the GCD;'
  write ( *, * ) '  I_TO_FAC converts an integer;'
  write ( *, * ) '  FAC_TO_I converts to an integer.'
  write ( *, * ) '  FAC_TO_RAT converts to a ratio.'

  i1 = 720
  i2 = 42

  call i_to_fac ( i1, nprime, npower1 )

  write ( *, * ) ' '
  write ( *, * ) '  Representation of I1 = ', i1
  write ( *, * ) ' '

  call fac_print ( nprime, npower1 )

  call i_to_fac ( i2, nprime, npower2 )

  write ( *, * ) ' '
  write ( *, * ) '  Representation of I2 = ', i2
  write ( *, * ) ' '

  call fac_print ( nprime, npower2 )

  call fac_lcm ( nprime, npower1, npower2, npower3 )

  write ( *, * ) ' '
  write ( *, * ) '  LCM of I1, I2:'
  write ( *, * ) ' '

  call fac_print ( nprime, npower3 )

  call fac_gcd ( nprime, npower1, npower2, npower3 )

  write ( *, * ) ' '
  write ( *, * ) '  GCD of I1, I2:'
  write ( *, * ) ' '

  call fac_print ( nprime, npower3 )

  call fac_mul ( nprime, npower1, npower2, npower3 )

  write ( *, * ) ' '
  write ( *, * ) '  Product of I1, I2:'
  write ( *, * ) ' '

  call fac_print ( nprime, npower3 )

  call fac_div ( nprime, npower2, npower1, npower3 )

  write ( *, * ) ' '
  write ( *, * ) '  Quotient of I2 / I1:'
  write ( *, * ) ' '

  call fac_print ( nprime, npower3 )

  call fac_to_rat ( nprime, npower3, itop, ibot )

  write ( *, * ) ' '
  write ( *, * ) 'Quotient as a rational: ', itop, ' / ', ibot

  return
end
subroutine test012
!
!*******************************************************************************
!
!! TEST012 tests GET_SEED.
!
  integer i
  integer iseed
  integer iseed_0
  integer iseed_1
  integer iseed_2
  integer iseed_3
  real uniform_01_sample
  real x

  write ( *, * ) ' '
  write ( *, * ) 'TEST012'
  write ( *, * ) '  GET_SEED gets a seed for the random number'
  write ( *, * ) '    generator.  These values are computed from'
  write ( *, * ) '    the time and date.  Values computed nearby'
  write ( *, * ) '    in time will be near to each other, and'
  write ( *, * ) '    should be passed through a random number'
  write ( *, * ) '    generator a few times before use.'
  write ( *, * ) ' '
  write ( *, * ) '     I         R(I)        R2(I)        R3(I)'
  write ( *, * ) ' '
  do i = 1, 10
    call get_seed ( iseed )
    iseed_0 = iseed
    x = uniform_01_sample ( iseed )
    iseed_1 = iseed
    x = uniform_01_sample ( iseed )
    iseed_2 = iseed
    x = uniform_01_sample ( iseed )
    iseed_3 = iseed
    write ( *, '(4i12)' ) iseed_0, iseed_1, iseed_2, iseed_3
  end do

  return
end
subroutine test013
!
!*******************************************************************************
!
!! TEST013 tests GRID1.
!
  integer, parameter :: ndim = 5
  integer, parameter :: nstep = 11
!
  integer j
  real x(ndim,nstep)
  real x1(ndim)
  real x2(ndim)
!
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST013'
  write ( *, * ) '  GRID1 computes a 1D grid between'
  write ( *, * ) '    two NDIM dimensional points X1 and X2.'
  write ( *, * ) ' '
  write ( *, * ) '  Here, we will use ', nstep, ' steps'
  write ( *, * ) '  going from '
  write ( *, '(5g12.4)' ) ( x1(j), j = 1, ndim )
  write ( *, * ) '  to'
  write ( *, '(5g12.4)' ) ( x2(j), j = 1, ndim )
  write ( *, * ) ' '
 
  call grid1 ( ndim, nstep, x, x1, x2 )
 
  call rmat_print ( ndim, ndim, nstep, x, '  The grid matrix:' )
 
  return
end
subroutine test014
!
!*******************************************************************************
!
!! TEST014 tests GRID1N.
!
  integer, parameter :: ndim = 5
  integer, parameter :: nstep = 11
!
  integer i
  integer j
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
 
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST014'
  write ( *, * ) '  GRID1N computes a 1D grid between'
  write ( *, * ) '    two NDIM dimensional points X1 and X2,'
  write ( *, * ) '    one point at a time.'
  write ( *, * ) ' '
  write ( *, * ) '  Here, we will use ', nstep, ' steps'
  write ( *, * ) '  going from '
  write ( *, '(5g12.4)' ) ( x1(j), j = 1, ndim )
  write ( *, * ) '  to'
  write ( *, '(5g12.4)' ) ( x2(j), j = 1, ndim )
  write ( *, * ) ' '
 
  do i = 1, nstep
    call grid1n ( i, ndim, nstep, x, x1, x2 )
    write ( *, '(i3,5g12.4)' ) i, ( x(j), j = 1, ndim )
  end do
 
  return
end
subroutine test015
!
!*******************************************************************************
!
!! TEST015 tests GRID2.
!
  integer, parameter :: ndim = 5
  integer, parameter :: nstep = 20
!
  integer i1
  integer i2
  integer j
  real x(ndim,nstep)
  real x1(ndim)
  real x2(ndim)
!
  i1 = 3
  i2 = 13
 
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST015'
  write ( *, * ) '  GRID2 computes a 1 D grid between'
  write ( *, * ) '    two NDIM dimensional points X1 and X2,'
  write ( *, * ) '    computing X1 and X2 at user specified times.'
  write ( *, * ) ' '
  write ( *, * ) '  Here, we will use ', nstep, ' steps,'
  write ( *, * ) '  and on step ', i1, ' we will compute'
  write ( *,'(5g12.4)' ) ( x1(j), j = 1, ndim )
  write ( *, * ) '  and on step ', i2, ' we will compute'
  write ( *,'(5g12.4)' ) ( x2(j), j = 1, ndim )
  write ( *, * ) ' '
 
  call grid2 ( i1, i2, ndim, nstep, x, x1, x2 )
 
  call rmat_print ( ndim, ndim, nstep, x, '  The grid matrix:' )

  return
end
subroutine test016
!
!*******************************************************************************
!
!! TEST016 tests GRID2N.
!
  integer, parameter :: ndim = 5
!
  integer i
  integer i1
  integer i2
  integer j
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
!
  i1 = 3
  i2 = 13
 
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST016'
  write ( *, * ) '  GRID2N computes points from a 1D grid'
  write ( *, * ) '    between two NDIM dimensional points'
  write ( *, * ) '    X1 and X2, one at a time, with X1 and X2'
  write ( *, * ) '    having user specified I coordinates.'
  write ( *, * ) ' '
  write ( *, * ) '  Here, on step ', i1, ' we would compute'
  write ( *, '(5g12.4)' ) ( x1(j), j = 1, ndim )
  write ( *, * ) '  and on step ', i2, ' we would compute'
  write ( *, '(5g12.4)' ) ( x2(j), j = 1, ndim )
  write ( *, * ) ' '
 
  do i = 1, 20
    call grid2n ( i, i1, i2, ndim, x, x1, x2 )
    write ( *, '(i3,5g12.4)' ) i, ( x(j), j = 1, ndim )
  end do
 
  return
end
subroutine test017
!
!*******************************************************************************
!
!! TEST017 tests GRID3.
!
  integer, parameter :: ndim = 5
  integer, parameter :: nstep1 = 3
  integer, parameter :: nstep2 = 6
!
  integer i
  integer j
  integer k
  real x(ndim,nstep1,nstep2)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  x3(1) = 1.0E+00
  x3(2) = 5.0E+00
  x3(3) = 0.0E+00
  x3(4) = 0.0E+00
  x3(5) = 3.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST017'
  write ( *, * ) '  GRID3 computes a 2D grid in the plane'
  write ( *, * ) '    containing the NDIM-dimensional'
  write ( *, * ) '    points X1, X2 and X3.'
  write ( *, * ) ' '
  write ( *, * ) '  Here, we will use ', nstep1, ' steps'
  write ( *, * ) '  going from '
  write ( *, '(5g12.4)' ) ( x1(k), k = 1, ndim )
  write ( *, * ) '  to'
  write ( *, '(5g12.4)' ) ( x2(k), k = 1, ndim )
  write ( *, * ) '  and ', nstep2,' steps going to '
  write ( *, '(5g12.4)' ) ( x3(k), k = 1, ndim )
  write ( *, * ) ' '
 
  call grid3 ( ndim, nstep1, nstep2, x, x1, x2, x3 )
 
  do i = 1, nstep1
    write ( *, * ) ' '
    do j = 1, nstep2
 
      write ( *, '(i3,i3,5g12.4)' ) i, j, ( x(k,i,j), k = 1, ndim )
 
    end do
  end do
 
  return
end
subroutine test018
!
!*******************************************************************************
!
!! TEST018 tests GRID3N.
!
  integer, parameter :: ndim = 5
  integer, parameter :: nstep1 = 3
  integer, parameter :: nstep2 = 6
!
  integer i
  integer j
  integer k
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  x3(1) = 1.0E+00
  x3(2) = 5.0E+00
  x3(3) = 0.0E+00
  x3(4) = 0.0E+00
  x3(5) = 3.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST018'
  write ( *, * ) '  GRID3D computes a point from a 2D'
  write ( *, * ) '    grid in the plane containing the '
  write ( *, * ) '    NDIM-dimensional points X1, X2 and X3.'
  write ( *, * ) ' '
  write ( *, * ) '  We use ', nstep1, ' steps from '
  write ( *, '(5g12.4)' ) ( x1(k), k = 1, ndim )
  write ( *, * ) '  to'
  write ( *, '(5g12.4)' ) ( x2(k), k = 1, ndim )
  write ( *, * ) '  and ', nstep2, ' steps going to '
  write ( *, '(5g12.4)' ) ( x3(k), k = 1, ndim )
  write ( *, * ) ' '
 
  do i = 1, nstep1
    write ( *, * ) ' '
    do j = 1, nstep2
 
      call grid3n ( i, j, ndim, nstep1, nstep2, x, x1, x2, x3 )
      write ( *, '(i3,i3,5g12.4)' ) i, j, ( x(k), k = 1, ndim )
 
    end do
  end do
 
  return
end
subroutine test019
!
!*******************************************************************************
!
!! TEST019 tests GRID4.
!
  integer, parameter :: ndim = 5
  integer, parameter :: nstep1 = 6
  integer, parameter :: nstep2 = 10
!
  integer i
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer k
  real x(ndim,nstep1,nstep2)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  i1 = 2
  i2 = 5
  j1 = 3
  j2 = 9
 
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  x3(1) = 1.0E+00
  x3(2) = 5.0E+00
  x3(3) = 0.0E+00
  x3(4) = 0.0E+00
  x3(5) = 3.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST019'
  write ( *, * ) '  GRID4 computes a 2D planar grid'
  write ( *, * ) '    containing the NDIM-dimensional'
  write ( *, * ) '    points X1, X2 and X3.'
  write ( *, * ) ' '
  write ( *, * ) '  We compute the points on the following steps:'
  write ( *, * ) ' '
  write ( *, * ) '  X1 on step ', i1, j1
  write ( *, * ) '  X2 on step ', i2, j1
  write ( *, * ) '  X3 on step ', i1, j2
  write ( *, * ) ' '
  write ( *, * ) '  We use ', nstep1, ' steps in the I direction'
  write ( *, * ) '  and ', nstep2, ' steps in the J direction.'
  write ( *, * ) ' '
  write ( *, * ) '  The points X1, X2 and X3 are:'
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) ( x1(k), k = 1, ndim )
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) ( x2(k), k = 1, ndim )
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) ( x3(k), k = 1, ndim )
  write ( *, * ) ' '
 
  call grid4 ( i1, i2, j1, j2, ndim, nstep1, nstep2, x, x1, x2, x3 )
 
  do i = 1, nstep1
    write ( *, * ) ' '
    do j = 1, nstep2
 
      write ( *, '(i3,i3,5g12.4)' ) i, j, ( x(k,i,j), k = 1, ndim )
 
    end do
  end do
 
  return
end
subroutine test020
!
!*******************************************************************************
!
!! TEST020 tests GRID4N.
!
  integer, parameter :: ndim = 5
  integer, parameter :: nstep1 = 6
  integer, parameter :: nstep2 = 10
!
  integer i
  integer i1
  integer i2
  integer j
  integer j1
  integer j2
  integer k
  real x(ndim)
  real x1(ndim)
  real x2(ndim)
  real x3(ndim)
!
  i1 = 2
  i2 = 5
  j1 = 3
  j2 = 9
 
  x1(1) = 1.0E+00
  x1(2) = 0.0E+00
  x1(3) = 20.0E+00
  x1(4) = -5.0E+00
  x1(5) = 1.0E+00
 
  x2(1) = 1.0E+00
  x2(2) = 10.0E+00
  x2(3) = 0.0E+00
  x2(4) = 5.0E+00
  x2(5) = 2.0E+00
 
  x3(1) = 1.0E+00
  x3(2) = 5.0E+00
  x3(3) = 0.0E+00
  x3(4) = 0.0E+00
  x3(5) = 3.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST020'
  write ( *, * ) '  GRID4N computes, one at a time, points' 
  write ( *, * ) '    on a 2D grid in the plane containing'
  write ( *, * ) '    the NDIM-dimensional points X1, X2 and X3.'
  write ( *, * ) ' '
  write ( *, * ) '  We wish to compute the points on the following'
  write ( *, * ) '  steps:'
  write ( *, * ) ' '
  write ( *, * ) '  X1 on step ', i1, j1
  write ( *, * ) '  X2 on step ', i2, j1
  write ( *, * ) '  X3 on step ', i1, j2
  write ( *, * ) ' '
  write ( *, * ) '  We use ', nstep1, ' steps in the I direction'
  write ( *, * ) '  and ', nstep2, ' steps in the J direction.'
  write ( *, * ) ' '
  write ( *, * ) '  The points X1, X2 and X3 are:'
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) ( x1(k),k = 1, ndim )
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) ( x2(k), k = 1, ndim )
  write ( *, * ) ' '
  write ( *, '(5g12.4)' ) ( x3(k), k = 1, ndim )
  write ( *, * ) ' '
 
  do i = 1, nstep1
    write ( *, * ) ' '
    do j = 1, nstep2
 
      call grid4n ( i, i1, i2, j, j1, j2, ndim, nstep1, nstep2, &
        x, x1, x2, x3 )

      write ( *, '(i3,i3,5g12.4)' ) i, j, ( x(k), k = 1, ndim )
 
    end do
  end do
 
  return
end
subroutine test021
!
!*******************************************************************************
!
!! TEST021 tests I_FACTOR.
!
  integer, parameter  :: maxfactor = 10
!
  integer factor(maxfactor)
  integer i
  integer n
  integer nfactor
  integer nleft
  integer power(maxfactor)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST021'
  write ( *, * ) '  I_FACTOR factors an integer,'

  n = 2**2 * 17 * 37 
  write ( *, * ) ' '
  write ( *, * ) '  The integer is ', n
  call i_factor ( n, maxfactor, nfactor, factor, power, nleft )
  write ( *, * ) ' '
  write ( *, * ) '  Prime representation:'
  write ( *, * ) ' '
  write ( *, * ) '  I, FACTOR(I), POWER(I)'
  write ( *, * ) ' '
  if ( abs ( nleft ) /= 1 ) then
    write ( *, * ) 0, nleft, ' (UNFACTORED PORTION)'
  end if

  do i = 1, nfactor
    write ( *, * ) i, factor(i), power(i)
  end do
 
  return
end
subroutine test022
!
!*******************************************************************************
!
!! TEST022 tests I_IS_PRIME.
!
  integer i
  logical i_is_prime
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST022'
  write ( *, * ) '  I_IS_PRIME reports whether an integer is prime.'

  write ( *, * ) ' '
  write ( *, * ) '  I     I_IS_PRIME(I)'
  write ( *, * ) ' '

  do i = -2, 25
    write ( *, * ) i, i_is_prime ( i )
  end do

  return
end
subroutine test023
!
!*******************************************************************************
!
!! TEST023 tests I_LOG_2.
!! TEST023 tests R_LOG_2.
!
  integer, parameter :: n = 20
!
  integer i
  integer i_log_2
  real r_log_2
  real temp
  real x(n)
!
  x(1) = 0.0E+00
  x(2) = 1.0E+00
  x(3) = 2.0E+00
  x(4) = 3.0E+00
  x(5) = 9.0E+00
  x(6) = 10.0E+00
  x(7) = 11.0E+00
  x(8) = 99.0E+00
  x(9) = 101.0E+00
  x(10) = -1.0E+00
  x(11) = -2.0E+00
  x(12) = -3.0E+00
  x(13) = -9.0E+00
  x(14) = 1.0E+00 / 2.0E+00
  x(15) = 1.0E+00 / 3.0E+00
  x(16) = 1.0E+00 / 4.0E+00
  x(17) = 1.0E+00 / 5.0E+00
  x(18) = 1.0E+00 / 99.0E+00
  x(19) = 1.0E+00 / 100.0E+00
  x(20) = 1.0E+00 / 101.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST023'
  write ( *, * ) '  I_LOG_2: whole part of log base 2.'
  write ( *, * ) '  R_LOG_2: log base 2,'
  write ( *, * ) ' '
  write ( *, * ) '  X, LOG(X)/LOG(2), R_LOG_2, I_LOG_2'
  write ( *, * ) ' '

  do i = 1, n
 
    if ( x(i) > 0.0E+00 ) then
      temp = log ( x(i) ) / log ( 2.0E+00 )
    else
      temp = 0.0E+00
    end if

    write ( *, '( 3g14.6, i6 )' ) x(i), temp, r_log_2 ( x(i) ), &
      i_log_2 ( x(i) )
 
  end do

  return
end
subroutine test024
!
!*******************************************************************************
!
!! TEST024 tests I_LOG_10.
!
  integer, parameter :: n = 20
!
  integer i
  integer i_log_10
  real temp
  real x(n)
!
  x(1) = 0.0E+00
  x(2) = 1.0E+00
  x(3) = 2.0E+00
  x(4) = 3.0E+00
  x(5) = 9.0E+00
  x(6) = 10.0E+00
  x(7) = 11.0E+00
  x(8) = 99.0E+00
  x(9) = 101.0E+00
  x(10) = -1.0E+00
  x(11) = -2.0E+00
  x(12) = -3.0E+00
  x(13) = -9.0E+00
  x(14) = 1.0E+00 / 2.0E+00
  x(15) = 1.0E+00 / 3.0E+00
  x(16) = 1.0E+00 / 4.0E+00
  x(17) = 1.0E+00 / 5.0E+00
  x(18) = 1.0E+00 / 99.0E+00
  x(19) = 1.0E+00 / 100.0E+00
  x(20) = 1.0E+00 / 101.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST024'
  write ( *, * ) '  I_LOG_10: whole part of log base 10,'
  write ( *, * ) ' '
  write ( *, * ) '  X, LOG10(|X|), I_LOG_10'
  write ( *, * ) ' '

  do i = 1, n
 
    if ( x(i) == 0.0E+00 ) then
      temp = 0.0E+00
    else
      temp = log10 ( abs ( x(i) ) )
    end if

    write ( *, '( 2g14.6, i6 )' ) x(i), temp, i_log_10 ( x(i) )
 
  end do

  return
end
subroutine test025
!
!*******************************************************************************
!
!! TEST025 tests I_MODDIV;
!! TEST025 tests I_MODP.
!
  integer, parameter :: ntest = 4
!
  integer i
  integer i_modp
  integer ndivid(ntest)
  integer nmult
  integer nrem
  integer number(ntest)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST025'
  write ( *, * ) '  I_MODDIV factors a number'
  write ( *, * ) '    into a multiple and a remainder.'
  write ( *, * ) ' '
  write ( *, * ) '    Number   Divisor  Multiple Remainder'
  write ( *, * ) ' '
 
  number(1) = 107
  ndivid(1) = 50
  number(2) = 107
  ndivid(2) = -50
  number(3) = -107
  ndivid(3) = 50
  number(4) = -107
  ndivid(4) = -50
 
  do i = 1, ntest
    call i_moddiv ( number(i), ndivid(i), nmult, nrem )
    write ( *, '(4i10)' ) number(i), ndivid(i), nmult, nrem
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Repeat using FORTRAN MOD:'
  write ( *, * ) ' '

  do i = 1, ntest
    nrem = mod ( number(i), ndivid(i) )
    nmult = number(i) / ndivid(i)
    write ( *, '(4i10)' ) number(i), ndivid(i), nmult, nrem
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Repeat using I_MODP:'
  write ( *, * ) ' '

  do i = 1, ntest
    nrem = i_modp ( number(i), ndivid(i) )
    nmult = ( number(i) - nrem ) / ndivid(i)
    write ( *, '(4i10)' ) number(i), ndivid(i), nmult, nrem
  end do
 
  return
end
subroutine test026
!
!*******************************************************************************
!
!! TEST026 tests I_ROUNDUP.
!
  integer i
  integer i_roundup
  integer ival
  real rval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST026'
  write ( *, * ) '  I_ROUNDUP rounds reals up.'
  write ( *, * ) ' '

  do i = -6, 6
    rval = real ( i ) / real ( 5.0 )
    ival = i_roundup ( rval )
    write ( *, '(g14.6,i8)' ) rval, ival
  end do

  return
end
subroutine test027
!
!*******************************************************************************
!
!! TEST027 tests I_GCD.
!! TEST027 tests I_LCM.
!
  integer, parameter :: ntest = 7
!
  integer i
  integer i_gcd
  integer itest
  integer ival(ntest)
  integer j
  integer jval(ntest)
  integer i_lcm
!
  ival(1) = 36
  ival(2) = 49
  ival(3) =  0
  ival(4) = 12
  ival(5) = 36
  ival(6) =  1
  ival(7) = 91

  jval(1) = 30
  jval(2) = -7
  jval(3) = 71
  jval(4) = 12
  jval(5) = 49
  jval(6) = 42
  jval(7) = 28

  write ( *, * ) ' '
  write ( *, * ) 'TEST027'
  write ( *, * ) '  I_GCD computes the greatest common factor,'
  write ( *, * ) '  I_LCM computes the least common multiple.'
  write ( *, * ) ' '
  write ( *, * ) '     I     J   I_GCD   I_LCM'
  write ( *, * ) ' '
 
  do itest = 1, ntest
    i = ival(itest)
    j = jval(itest)
    write ( *, '(4i6)') i, j, i_gcd(i,j), i_lcm(i,j)
  end do
 
  return
end
subroutine test028
!
!*******************************************************************************
!
!! TEST028 tests I_JACOBI_SYMBOL.
!
  integer, parameter :: ntest = 4
!
  integer i
  integer l
  integer p
  integer ptest(ntest)
  integer q
!
  ptest(1) = 3
  ptest(2) = 9
  ptest(3) = 10
  ptest(4) = 12

  write ( *, * ) ' '
  write ( *, * ) 'TEST028'
  write ( *, * ) '  I_JACOBI_SYMBOL computes the Jacobi symbol'
  write ( *, * ) '    (Q/P), which records if Q is a quadratic '
  write ( *, * ) '    residue modulo the number P.'

  do i = 1, ntest
    p = ptest(i)
    write ( *, * ) ' '
    write ( *, * ) 'Jacobi Symbols for P = ', p
    write ( *, * ) ' '
    do q = 0, p
      call i_jacobi_symbol ( q, p, l )
      write ( *, '(3i8)' ) p, q, l 
    end do
  end do

  return
end
subroutine test029
!
!*******************************************************************************
!
!! TEST029 tests I_LEGENDRE_SYMBOL.
!
  integer, parameter :: ntest = 4
!
  integer i
  integer l
  integer p
  integer ptest(ntest)
  integer q
!
  ptest(1) = 7
  ptest(2) = 11
  ptest(3) = 13
  ptest(4) = 17

  write ( *, * ) ' '
  write ( *, * ) 'TEST029'
  write ( *, * ) '  I_LEGENDRE_SYMBOL computes the Legendre'
  write ( *, * ) '    symbol (Q/P) which records whether Q is '
  write ( *, * ) '    a quadratic residue modulo the prime P.'

  do i = 1, ntest
    p = ptest(i)
    write ( *, * ) ' '
    write ( *, * ) 'Legendre Symbols for P = ', p
    write ( *, * ) ' '
    do q = 0, p
      call i_legendre_symbol ( q, p, l )
      write ( *, '(3i8)' ) p, q, l 
    end do
  end do

  return
end
subroutine test030
!
!*******************************************************************************
!
!! TEST030 tests I_MANT.
!
  integer is
  integer j
  integer k
  integer l
  real x
!
  x = -314.159E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST030'
  write ( *, * ) '  I_MANT decomposes an integer,'
  write ( *, * ) ' '
  write ( *, * ) '  Number to be decomposed is X = ', x

  call i_mant ( x, is, j, k, l )

  write ( *, * ) ' '
  write ( *, * ) '  I_MANT: X = ', is, ' * (', j, '/', k, ') * 2**', l
 
  return
end
subroutine test031
!
!*******************************************************************************
!
!! TEST031 tests I_MEMORY.
!
  integer inc
  integer ival
  integer ivan
  integer jack
  integer kyle
  integer lou
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST031'
  write ( *, * ) '  I_MEMORY'
  write ( *, * ) ' '

  ival = 0
  call i_memory ( 'INIT', ' ', ival )

  ivan = 1
  jack = 2
  kyle = 3

  call i_memory ( 'NAME', 'IVAN', ivan )
  call i_memory ( 'NAME', 'JACK', jack )
  call i_memory ( 'NAME', 'KYLE', kyle )

  call i_memory ( 'PRINT', '*', ival )

  write ( *, * ) ' '
  inc = 100
  call i_memory ( 'PRINT', 'IVAN', ivan )
  write ( *, * ) 'INC IVAN 100:'
  call i_memory ( 'INC', 'IVAN', inc )
  call i_memory ( 'PRINT', 'IVAN', ivan )

  lou = -33
  call i_memory ( 'NAME', 'LOU', lou )

  write ( *, * ) ' '
  write ( *, * ) 'PUSH IVAN 17'
  ivan = 17
  call i_memory ( 'PUSH', 'IVAN', ivan )
  call i_memory ( 'PRINT', 'IVAN', ivan )
 
  write ( *, * ) 'INC IVAN 100'
  call i_memory ( 'INC', 'IVAN', inc )
  call i_memory ( 'PRINT', 'IVAN', ivan )
 
  write ( *, * ) 'POP IVAN = ',ivan
  call i_memory ( 'POP', 'IVAN', ivan )
  call i_memory ( 'PRINT', 'IVAN', ivan )

  write ( *, * ) ' '
  write ( *, * ) 'SET JACK 99:'
 
  jack = 99
  call i_memory ( 'SET', 'JACK', jack )
  call i_memory ( 'PRINT', 'JACK', jack )
 
  write ( *, * ) ' '
  write ( *, * ) 'INC JACK 100:'
 
  call i_memory ( 'INC', 'JACK', inc )
  call i_memory ( 'PRINT', 'JACK', jack )

  write ( *, * ) ' '
  write ( *, * ) 'PRINT *'
  call i_memory ( 'PRINT', '*', ival )
 
  return
end
subroutine test032
!
!*******************************************************************************
!
!! TEST032 tests I_MOEBIUS;
!! TEST032 tests I_OMEGA;
!! TEST032 tests I_PHI;
!! TEST032 tests I_SIGMA;
!! TEST032 tests I_TAU.
!
  integer mu
  integer n
  integer omega
  integer phin
  integer sigman
  integer taun
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST032'
  write ( *, * ) '  I_MOEBIUS computes the Moebius function;'
  write ( *, * ) '  I_OMEGA counts the distinct prime divisors of N;'
  write ( *, * ) '  I_PHI computes the number of values less than N'
  write ( *, * ) '    and relatively prime to N;'
  write ( *, * ) '  I_TAU computes the number of factors;'
  write ( *, * ) '  I_SIGMA computes the sum of factors.'
  write ( *, * ) ' '
  write ( *, * ) '  N  PHI(N), TAU(N), SIGMA(N), MU(N)  OMEGA(N)'
  write ( *, * ) ' '
  do n = 1, 20
    call i_moebius ( n, mu )
    call i_omega ( n, omega )
    call i_phi ( n, phin )
    call i_tau ( n, taun )
    call i_sigma ( n, sigman )
    write ( *, '(6i5)' ) n, phin, taun, sigman, mu, omega
  end do

  return
end
subroutine test033
!
!*******************************************************************************
!
!! TEST033 tests I_SIGN.
!
  integer, parameter :: test_num = 5
!
  integer test_i
  integer i_sign
  integer x
  integer x_test(test_num)
!
  x_test(1) = - 10
  x_test(2) = - 7
  x_test(3) =   0
  x_test(4) = + 5
  x_test(5) = + 9

  write ( *, * ) ' '
  write ( *, * ) 'TEST033'
  write ( *, * ) '  I_SIGN returns the sign of a number.'
  write ( *, * ) ' '
 
  do test_i = 1, test_num
    x = x_test(test_i)
    write ( *, '(2i6)' ) x, i_sign ( x )
  end do
 
  return
end
subroutine test888
!
!*******************************************************************************
!
!! TEST888 tests R_TO_R_DISCRETE.
!
  integer i
  integer ndx
  real r
  real rd
  real rhi
  real rlo
!
  rlo = 1.0E+00
  rhi = 10.0E+00
  ndx = 19

  write ( *, * ) ' '
  write ( *, * ) 'TEST888'
  write ( *, * ) '  R_TO_R_DISCRETE maps real numbers to a discrete set'
  write ( *, * ) '  of equally spaced real numbers in an interval.'
  write ( *, * ) ' '
  write ( *, * ) '  Number of discrete values = ', ndx
  write ( *, * ) '  Real interval: ', rlo, rhi
  write ( *, * ) ' '
  write ( *, * ) '  R   RD'
  write ( *, * ) ' '

  do i = 1, 15
    call r_random ( rlo-2.0, rhi+2.0, r )
    call r_to_r_discrete ( r, rlo, rhi, ndx, rd )
    write ( *, '(g14.6,g14.6)' ) r, rd
  end do

  return
end
subroutine test034
!
!*******************************************************************************
!
!! TEST034 tests IINT_TO_RINT;
!! TEST034 tests RINT_TO_IINT;
!
  integer i
  integer ihi
  integer ilo
  integer ir
  real r
  real r2
  real rhi
  real rlo
!
  rlo = 100.0E+00
  rhi = 200.0E+00
  ilo = 1
  ihi = 11

  write ( *, * ) ' '
  write ( *, * ) 'TEST034'
  write ( *, * ) '  For data in an interval,'
  write ( *, * ) '  IINT_TO_RINT converts an integer to a real;'
  write ( *, * ) '  RINT_TO_IINT converts a real to an integer.'
  write ( *, * ) ' '
  write ( *, * ) '  Integer interval: ', ilo, ihi
  write ( *, * ) '  Real interval: ', rlo, rhi
  write ( *, * ) ' '
  write ( *, * ) '  R   I(R)  R(I(R))'
  write ( *, * ) ' '

  do i = 1, 10
    call r_random ( rlo-15.0, rhi+15.0, r )
    call rint_to_iint ( rlo, rhi, r, ilo, ihi, ir )
    call iint_to_rint ( ilo, ihi, ir, rlo, rhi, r2 )
    write ( *, '(g14.6,i6,g14.6)' ) r, ir, r2
  end do

  return
end
subroutine test0345
!
!*******************************************************************************
!
!! TEST0345 tests I_WRAP.
!
  integer i
  integer i_wrap
  integer ihi
  integer ilo
!
  ilo = 4
  ihi = 8

  write ( *, * ) ' '
  write ( *, * ) 'TEST0345'
  write ( *, * ) '  I_WRAP forces an integer to lie within given limits.'
  write ( *, * ) ' '
  write ( *, * ) '  ILO = ', ilo, ' IHI = ', ihi
  write ( *, * ) ' '
  write ( *, * ) '     I  I_WRAP(I)'
  write ( *, * ) ' '

  do i = -10, 20
    write ( *, '(2i6)' ) i, i_wrap ( i, ilo, ihi )
  end do

  return
end
subroutine test035
!
!*******************************************************************************
!
!! TEST035 tests ICOL_SORT_A.
!! TEST035 tests ICOL_SORT_D.
!
  integer, parameter :: m = 5
  integer, parameter :: n = 4
  integer, parameter :: lda = m
!
  integer a(lda,n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST035'
  write ( *, * ) '  ICOL_SORT_A ascending sorts an integer array'
  write ( *, * ) '    as a table of columns.'
  write ( *, * ) '  ICOL_SORT_D descending sorts an integer array'
  write ( *, * ) '    as a table of columns.'

  call imat_random ( 1, 10, lda, m, n, a )
 
  call imat_print ( lda, m, n, a, '  The original matrix:' )
 
  call icol_sort_a ( lda, m, n, a )
 
  call imat_print ( lda, m, n, a, '  Ascending sorted:' )

  call icol_sort_d ( lda, m, n, a )
 
  call imat_print ( lda, m, n, a, '  Descending sorted:' )

  return
end
subroutine test036
!
!*******************************************************************************
!
!! TEST036 tests IMAT_ELIM.
!! TEST036 tests IMAT_RED.
!
  integer, parameter :: lda = 5
!
  integer a(lda,lda)
  integer i
  integer icol(lda)
  integer ifact
  integer imat
  integer irow(lda)
  integer j
  integer k
  integer m
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST036'
  write ( *, * ) '  IMAT_RED divides common factors in a matrix;'
  write ( *, * ) '  IMAT_ELIM does exact Gauss elimination.'
  write ( *, * ) ' '
 
  do imat = 1, 3
 
    if ( imat == 1 ) then
 
      m = 5
      n = 5
 
      k = 0
      do i = 1, m
        do j = 1, n
          k = k+1
          a(i,j) = k
        end do
      end do
 
    else if ( imat == 2 ) then
 
      m = 5
      n = 5
 
      ifact = 8 * 7 * 6 * 5 * 4 * 3 * 2
 
      do i = 1, m
        do j = 1, n
          a(i,j) = ifact / ( i + j - 1 )
        end do
      end do
 
    else if ( imat == 3 ) then
 
      m = 4
      n = 5
 
      do i = 1, m
        do j = 1, n
          a(i,j) = i * j
        end do
      end do

    end if
 
    call imat_print ( lda, m, n, a, '  The original matrix:' )
 
    call imat_red ( lda, m, n, a, irow, icol )
 
    write ( *, * ) ' '
    write ( *, * ) '  The matrix, as returned by IMAT_RED:'
    write ( *, * ) ' '
    do i = 1, m
      write ( *, '(6i8)' )  ( a(i,j), j = 1, n ), irow(i)
    end do
    write ( *, '(5i8)' )  ( icol(j), j = 1, n )
 
    call imat_elim ( lda, n, n, a )
 
    call imat_print ( lda, m, n, a, '  The matrix returned by IMAT_ELIM:' )
 
  end do
 
  return
end
subroutine test0365
!
!*******************************************************************************
!
!! TEST0365 tests IMAT_L1_INVERSE.
!
  integer, parameter :: n = 6
  integer, parameter :: lda = n
!
  integer a(lda,n)
  integer b(lda,n)
  integer c(lda,n)
  integer i
  integer j
  integer k
!
  data ( ( a(i,j), j = 1, n ), i = 1, n ) / &
     1, 0, 0, 0, 0, 0, &
     2, 1, 0, 0, 0, 0, &
     0, 0, 1, 0, 0, 0, &
     5, 0, 3, 1, 0, 0, &
     0, 0, 0, 0, 1, 0, &
    75, 0, 0, 6, 4, 1 /
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0365'
  write ( *, * ) '  IMAT_L1_INVERSE inverts a unit lower triangular matrix.'
  write ( *, * ) ' '

  call imat_print ( lda, n, n, a, '  The original matrix:' )
 
  call imat_l1_inverse ( lda, n, a, b )
 
  call imat_print ( lda, n, n, b, '  The inverse matrix:' )
 
  do i = 1, n
    do j = 1, n
      c(i,j) = 0.0
      do k = 1, n
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do

  call imat_print ( lda, n, n, c, '  The product:' )

  return
end
subroutine test0366
!
!*******************************************************************************
!
!! TEST0366 tests IMAT_IMAX.
!! TEST0366 tests IMAT_IMIN.
!! TEST0366 tests IMAT_MAX.
!! TEST0366 tests IMAT_MIN.
!
  integer, parameter :: m = 5
  integer, parameter :: lda = m
  integer, parameter :: n = 7
!
  integer a(lda,n)
  integer i
  integer j
  integer imat_max
  integer imat_min
  integer temp1
  integer temp2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0366'
  write ( *, * ) '  IMAT_IMAX locates the maximum;'
  write ( *, * ) '  IMAT_IMIN locates the maximum;'
  write ( *, * ) '  IMAT_MAX computes the maximum;'
  write ( *, * ) '  IMAT_MIN computes the minimum;'
  write ( *, * ) ' '
 
  call imat_random ( 0, 10, lda, m, n, a )
 
  call imat_print ( lda, m, n, a, '  Random array:' )
 
  temp1 = imat_min ( lda, m, n, a )
  temp2 = imat_max ( lda, m, n, a )

  write ( *, * ) ' '
  write ( *, * ) '  Minimum, Maximum =             ', temp1, temp2
  call imat_imax ( lda, m, n, a, i, j )
  write ( *, * ) '  Maximum I,J indices            ', i, j
  call imat_imin ( lda, m, n, a, i, j )
  write ( *, * ) '  Minimum I,J indices            ', i, j

  return
end
subroutine test037
!
!*******************************************************************************
!
!! TEST037 tests IMAT_PERM_RANDOM.
!
  integer, parameter :: lda = 5
  integer, parameter :: n = 5
!
  integer a(lda,n)
  integer i
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST037'
  write ( *, * ) '  IMAT_PERM_RANDOM applies a random permutation'
  write ( *, * ) '    to a square integer matrix.'
  write ( *, * ) ' '
 
  do i = 1, n
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call imat_print ( lda, n, n, a, '  The original matrix:' )

  call imat_perm_random ( lda, n, a )

  call imat_print ( lda, n, n, a, '  The permuted matrix:' )

  return
end
subroutine test038
!
!*******************************************************************************
!
!! TEST038 tests IMAT_PERM2_RANDOM.
!! TEST038 tests IROW_SORT_D;
!! TEST038 tests IROW_SORT2_D;
!
  integer, parameter :: m = 6
  integer, parameter :: lda = m
  integer, parameter :: n = 4
!
  integer a(lda,n)
  integer i
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST038'
  write ( *, * ) '  For a rectangular integer matrix:'
  write ( *, * ) '  IMAT_PERM2_RANDOM applies independent'
  write ( *, * ) '    random permutations to rows and columns;'
  write ( *, * ) '  IROW_SORT_D sorts the rows;'
  write ( *, * ) '  IROW_SORT2_D sorts the elements of the rows.'
 
  do i = 1, m
    do j = 1, n
      a(i,j) = 10 * i + j
    end do
  end do

  call imat_print ( lda, m, n, a, '  The original matrix:' )

  write ( *, * ) ' '

  call imat_perm2_random ( lda, m, n, a )

  call imat_print ( lda, m, n, a, '  The permuted matrix:' )

  call irow_sort_d ( lda, m, n, a )

  call imat_print ( lda, m, n, a, '  The row-sorted matrix:' )

  call irow_sort2_d ( lda, m, n, a )

  call imat_print ( lda, m, n, a, '  The element-sorted row-sorted matrix:' )

  return
end
subroutine test0385
!
!*******************************************************************************
!
!! TEST0385 tests IMAT_U1_INVERSE.
!
  integer, parameter :: n = 6
  integer, parameter :: lda = n
!
  integer a(lda,n)
  integer b(lda,n)
  integer c(lda,n)
  integer i
  integer j
  integer k
!
  data ( ( a(i,j), j = 1, n ), i = 1, n ) / &
    1, 2, 0, 5, 0, 75, &
    0, 1, 0, 0, 0,  0, &
    0, 0, 1, 3, 0,  0, &
    0, 0, 0, 1, 0,  6, &
    0, 0, 0, 0, 1,  4, &
    0, 0, 0, 0, 0,  1 /
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0385'
  write ( *, * ) '  IMAT_U1_INVERSE inverts a unit upper triangular matrix.'
  write ( *, * ) ' '

  call imat_print ( lda, n, n, a, '  The original matrix:' )
 
  call imat_u1_inverse ( lda, n, a, b )
 
  call imat_print ( lda, n, n, b, '  The inverse matrix:' )

  do i = 1, n
    do j = 1, n
      c(i,j) = 0.0E+00
      do k = 1, n
        c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
    end do
  end do

  call imat_print ( lda, n, n, c, '  The product:' )

  return
end
subroutine test039
!
!*******************************************************************************
!
!! TEST039 tests IROW_MAX;
!! TEST039 tests IROW_MEAN;
!! TEST039 tests IROW_MIN;
!! TEST039 tests IROW_SUM;
!! TEST039 tests IROW_SWAP;
!! TEST039 tests IROW_VARIANCE.
!
  integer, parameter :: m = 3
  integer, parameter :: n = 4
  integer, parameter :: lda = m
!
  integer a(m,n)
  integer amax(m)
  integer amin(m)
  integer i
  integer irow1
  integer irow2
  integer imax(m)
  integer imin(m)
  integer j
  integer k
  real mean(m)
  integer sum(m)
  real variance(m)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST039'
  write ( *, * ) '  IROW_MAX computes row maximums;'
  write ( *, * ) '  IROW_MEAN computes row means;'
  write ( *, * ) '  IROW_MIN computes row minimums;'
  write ( *, * ) '  IROW_SUM computes row sums;'
  write ( *, * ) '  IROW_SWAP swaps two rows;'
  write ( *, * ) '  IROW_VARIANCE computes row variances;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = k
    end do
  end do

  call imat_print ( lda, m, n, a, '  The matrix:' )

  call irow_sum ( lda, m, n, a, sum )

  call irow_max ( lda, m, n, a, imax, amax )

  call irow_min ( lda, m, n, a, imin, amin )

  call irow_mean ( lda, m, n, a, mean )

  call irow_variance ( lda, m, n, a, variance )

  write ( *, * ) ' '
  write ( *, * ) 'Maximum, minimum, sum, mean, variance:'
  write ( *, * ) ' '
  do i = 1, m
    write ( *, '(i3,3x,3i6,2f10.4)' ) &
      i, amax(i), amin(i), sum(i), mean(i), variance(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Swap rows 1 and 3:'
  write ( *, * ) ' '

  irow1 = 1
  irow2 = 3
  call irow_swap ( lda, m, n, a, irow1, irow2 )

  call imat_print ( lda, m, n, a, '  The new matrix:' )

  return
end
subroutine test040
!
!*******************************************************************************
!
!! TEST040 tests ISBN_CHECK.
!
  integer, parameter :: ntest = 8
!
  integer check
  integer i
  character ( len = 20 ) isbn_test(ntest)
!
  isbn_test(1) = '0-8493-9640-9'
  isbn_test(2) = '0-201-54275-7'
  isbn_test(3) = '0-521-35796-9'
  isbn_test(4) = '0-07-034025-0'

  isbn_test(5) = '0-7493-9640-9'
  isbn_test(6) = '0-201-54275-X'
  isbn_test(7) = '0-521-X5796-9'
  isbn_test(8) = '0-37-034025-0'

  write ( *, * ) ' '
  write ( *, * ) 'TEST040'
  write ( *, * ) '  ISBN_CHECK checks ISBN''s.'
  write ( *, * ) ' '
  write ( *, * ) '  A correct ISBN has a checksum of 0.'
  write ( *, * ) ' '

  do i = 1, ntest

    call isbn_check ( isbn_test(i), check )

    write ( *, '(a20,5x,i2)' ) isbn_test(i), check

  end do

  return
end
subroutine test041
!
!*******************************************************************************
!
!! TEST041 tests ISBN_FILL.
!
  integer, parameter :: ntest = 5
!
  integer check
  integer i
  character ( len = 20 ) isbn
  character ( len = 20 ) isbn_test(ntest)
!
  isbn_test(1) = '0-?493-9640-9'
  isbn_test(2) = '0-201-5427?-7'
  isbn_test(3) = '0-521-35796-?'
  isbn_test(4) = '?-07-034025-0'
  isbn_test(5) = '0-07-05?489-2'

  write ( *, * ) ' '
  write ( *, * ) 'TEST041'
  write ( *, * ) '  ISBN_FILL can fill in a single missing digit'
  write ( *, * ) '  in an ISBN.'
  write ( *, * ) ' '

  do i = 1, ntest

    isbn = isbn_test(i) 

    call isbn_fill ( isbn )

    call isbn_check ( isbn, check )

    write ( *, '(a20,5x,a20,5x,i2)' ) isbn_test(i), isbn, check

  end do

  return
end
subroutine test042
!
!*******************************************************************************
!
!! TEST042 tests IVEC_AMAX;
!! TEST042 tests IVEC_AMIN;
!! TEST042 tests IVEC_AMINZ;
!! TEST042 tests IVEC_CUM;
!! TEST042 tests IVEC_IMAX;
!! TEST042 tests IVEC_IMAX_LAST;
!! TEST042 tests IVEC_INDEX;
!! TEST042 tests IVEC_MAX;
!! TEST042 tests IVEC_MEDIAN;
!! TEST042 tests IVEC_MEAN;
!! TEST042 tests IVEC_MIN;
!! TEST042 tests IVEC_NONZERO;
!! TEST042 tests IVEC_VARIANCE.
!
  integer, parameter :: n = 10
!
  integer a(n)
  integer a_cum(n+1)
  integer aval
  integer ival
  integer ivec_imax_last
  integer ivec_index
  integer ivec_nonzero
  integer j
  real mean
  integer median
  integer nonzero
  real variance
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST042'
  write ( *, * ) '  For an integer vector:'
  write ( *, * ) '  IVEC_AMAX:      maximum absolute entry;'
  write ( *, * ) '  IVEC_AMIN:      minimum absolute entry;'
  write ( *, * ) '  IVEC_AMINZ:     minimum nonzero absolute entry;'
  write ( *, * ) '  IVEC_CUM:       cumulative sum;'
  write ( *, * ) '  IVEC_IMAX:      a maximal index;'
  write ( *, * ) '  IVEC_IMAX_LAST: last maximal index;'
  write ( *, * ) '  IVEC_INDEX:     first index of given value;'
  write ( *, * ) '  IVEC_MAX:       maximum entry;'
  write ( *, * ) '  IVEC_MEAN:      mean value;'
  write ( *, * ) '  IVEC_MEDIAN:    median value;'
  write ( *, * ) '  IVEC_MIN:       minimum entry;'
  write ( *, * ) '  IVEC_NONZERO:   number of nonzeroes;'
  write ( *, * ) '  IVEC_VARIANCE:  variance.'
  write ( *, * ) ' '
 
  call ivec_random ( -n, n, n, a )
 
  call ivec_print ( n, a, '  Input vector:' )

  write ( *, * ) ' '

  call ivec_max ( n, a, aval )
  write ( *, * ) '  Maximum:                  ', aval
  call ivec_imax ( n, a, ival )
  write ( *, * ) '  Maximum index:            ', ival
  ival = ivec_imax_last ( n, a )
  write ( *, * ) '  Last maximum index:       ', ival

  call ivec_min ( n, a, aval )
  call ivec_imin ( n, a, ival )

  write ( *, * ) '  Minimum:                  ', aval
  write ( *, * ) '  Minimum index:            ', ival

  call ivec_amax ( n, a, aval )
  call ivec_iamax ( n, a, ival )

  write ( *, * ) '  Maximum abs:              ', aval
  write ( *, * ) '  Maximum abs index:        ', ival

  call ivec_amin ( n, a, aval )
  call ivec_iamin ( n, a, ival )

  write ( *, * ) '  Minimum abs:              ', aval
  write ( *, * ) '  Minimum abs index:        ', ival

  call ivec_aminz ( n, a, aval )
  call ivec_iaminz ( n, a, ival )

  write ( *, * ) '  Minimum abs nonzero:      ', aval
  write ( *, * ) '  Minimum abs nonzero index:', ival

  write ( *, * ) ' '
  j = ivec_index ( n, a, aval )
  write ( *, * ) '  Index of first occurrence of ', aval, ' is ', j

  aval = aval + 1
  j = ivec_index ( n, a, aval )
  write ( *, * ) '  Index of first occurrence of ', aval, ' is ', j

  call ivec_cum ( n, a, a_cum )
  call ivec_mean ( n, a, mean )
  call ivec_median ( n, a, median )
  nonzero = ivec_nonzero ( n, a )
  call ivec_variance ( n, a, variance )

  write ( *, * ) ' '
  write ( *, * ) '  Sum of entries:           ', sum ( a )
  write ( *, * ) '  Mean:                     ', mean
  write ( *, * ) '  Median:                   ', median
  write ( *, * ) '  Number of nonzeroes :     ', nonzero
  write ( *, * ) '  Variance:                 ', variance

  call ivec_print ( n+1, a_cum, '  Cumulative sums:' )
 
  return
end
subroutine test043
!
!*******************************************************************************
!
!! TEST043 tests IVEC_BRACKET.
!! TEST043 tests IVEC_INSERT.
!
  integer, parameter :: n_max = 20
  integer, parameter :: ntest = 6
!
  integer a(n_max)
  integer atest(ntest)
  integer aval
  integer i
  integer itest
  integer left
  integer n
  integer right
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST043'
  write ( *, * ) '  IVEC_BRACKET finds a pair of entries in a'
  write ( *, * ) '    sorted integer array which bracket a value.'
  write ( *, * ) '  IVEC_INSERT inserts a value into a vector.'
  write ( *, * ) ' '
  write ( *, * ) '  We use these two routines to bracket a value,'
  write ( *, * ) '  and then insert it.'

  atest(1) = -10
  atest(2) = 2
  atest(3) = 9
  atest(4) = 10
  atest(5) = 20
  atest(6) = 24

  n = 10
  do i = 1, n
    a(i) = 2 * i
  end do
  a(6) = a(5)

  call ivec_print ( n, a, '  Sorted array:' )

  do itest = 1, ntest

    aval = atest(itest)

    write ( *, * ) ' '
    write ( *, * ) 'Search for AVAL = ', aval
    call ivec_bracket ( n, a, aval, left, right )
    write ( *, * ) 'Left = ', left
    write ( *, * ) 'Right = ', right

    if ( left >= 1 ) then 
      write ( *, * ) 'A(LEFT)=', a(left)
    end if

    if ( right >= 1 ) then
      write ( *, * ) 'A(RIGHT) = ', a(right)
    end if
!
!  Insert the value.
!
    if ( left == -1 ) then
      left = 0
    end if

    if ( left == right ) then

      write ( *, * ) ' '
      write ( *, * ) 'No insertion necessary.'

    else

      call ivec_insert ( n, a, left+1, aval )

      n = n + 1

      call ivec_print ( n, a, '  Sorted, augmented array:' )

    end if

  end do

  return
end
subroutine test0435
!
!*******************************************************************************
!
!! TEST0435 tests IVEC_DESCENDS.
!
  integer, parameter :: n = 4
  integer, parameter :: ntest = 6
!
  integer itest
  logical ivec_descends
  integer x(n,ntest)
!
  x(1,1) = 1
  x(2,1) = 3
  x(3,1) = 2
  x(4,1) = 4

  x(1,2) = 2
  x(2,2) = 2
  x(3,2) = 2
  x(4,2) = 2

  x(1,3) = 1
  x(2,3) = 2
  x(3,3) = 2
  x(4,3) = 4

  x(1,4) = 1
  x(2,4) = 2
  x(3,4) = 3
  x(4,4) = 4

  x(1,5) = 4
  x(2,5) = 4
  x(3,5) = 3
  x(4,5) = 1

  x(1,6) = 9
  x(2,6) = 7
  x(3,6) = 3
  x(4,6) = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST0435'
  write ( *, * ) '  IVEC_DESCENDS determines if an integer vector descends.'
  write ( *, * ) ' '

  do itest = 1, ntest

    call ivec_print ( n, x(1,itest), '  Test vector:' )

    write ( *, * ) '  IVEC_DESCENDS = ', ivec_descends ( n, x(1,itest) )

  end do
   
  return
end
subroutine test044
!
!*******************************************************************************
!
!! TEST044 tests IVEC_FRAC;
!
  integer, parameter :: n = 10
!
  integer a(n)
  integer iafrac
  integer k
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST044'
  write ( *, * ) '  IVEC_FRAC: K-th smallest integer vector entry.'
  write ( *, * ) ' '

  call ivec_random ( 1, 2*n, n, a )

  call ivec_print ( n, a, '  The array to search:' )

  write ( *, * ) ' '
  write ( *, * ) '  Fractile    Value'
  write ( *, * ) ' '

  do k = 1, n, n/2

    call ivec_frac ( n, a, k, iafrac )

    write ( *, '(2i8)' ) k, iafrac

  end do

  return
end
subroutine test0445
!
!*******************************************************************************
!
!! TEST0445 tests IVEC_HEAP_A.
!! TEST0445 tests IVEC_HEAP_D.
!
  integer, parameter :: n = 10
!
  integer a(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0445'
  write ( *, * ) '  For an integer vector,'
  write ( *, * ) '  IVEC_HEAP_A puts into ascending heap form.'
  write ( *, * ) '  IVEC_HEAP_D puts into descending heap form.'
  write ( *, * ) ' '

  call ivec_random ( 0, n, n, a )

  call ivec_print ( n, a, '  Unsorted array:' )

  call ivec_heap_a ( n, a )
 
  call ivec_print ( n, a, '  Ascending heap form:' )
 
  call ivec_heap_d ( n, a )
 
  call ivec_print ( n, a, '  Descending heap form:' )

  return
end
subroutine test045
!
!*******************************************************************************
!
!! TEST045 tests IVEC_HEAP_D_EXTRACT.
!! TEST045 tests IVEC_HEAP_D_INSERT.
!! TEST045 tests IVEC_HEAP_D_MAX
!
  integer, parameter :: n_max = 10
!
  integer a(n_max)
  integer i
  integer n
  integer val
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST045'
  write ( *, * ) '  For a descending heap of integers,'
  write ( *, * ) '  IVEC_HEAP_D_INSERT inserts a value into the heap.'
  write ( *, * ) '  IVEC_HEAP_D_EXTRACT extracts the maximum value;'
  write ( *, * ) '  IVEC_HEAP_D_MAX reports the maximum value.'
  write ( *, * ) ' '
  write ( *, * ) '  These 3 operations are enough to model a priority queue.'

  n = 0

  do i = 1, n_max

    call i_random ( 0, 10, val )

    call ivec_heap_d_insert ( n, a, val )

    write ( *, * ) ' '
    write ( *, * ) '  Inserting value          ', val

    call ivec_heap_d_max ( n, a, val )

    write ( *, * ) '  Current maximum value is ', val

  end do

  call ivec_print ( n, a, '  Current heap as a vector:' )

  write ( *, * ) ' '
  write ( *, * ) '  Now extract the maximum several times.'
  write ( *, * ) ' '

  do i = 1, 5
    call ivec_heap_d_extract ( n, a, val )
    write ( *, * ) 'Extracting maximum element = ', val
  end do

  call ivec_print ( n, a, '  Current heap as a vector:' )

  return
end
subroutine test0451
!
!*******************************************************************************
!
!! TEST0451 tests IVEC_INDEX_SEARCH.
!! TEST0451 tests IVEC_INDEX_INSERT_UNIQUE.
!
  integer, parameter :: max_n = 20
!
  integer equal
  integer i
  integer indx(max_n)
  integer j
  integer less
  integer more
  integer n
  integer x(max_n)
  integer xval
!
  n = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST0451'
  write ( *, * ) '  IVEC_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, * ) '    index sorted array.'
  write ( *, * ) '  IVEC_INDEX_SEARCH searches for an entry with a given value.'
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values:'
  write ( *, * ) ' '
  do i = 1, 20
    call i_random ( 0, 20, xval )
    call ivec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Results of search for given XVAL:'
  write ( *, * ) ' '
  write ( *, '(a)' ) 'XVAL  Less Equal More'
  write ( *, * ) ' '

  do xval = 0, 20
    call ivec_index_search ( n, x, indx, xval, less, equal, more )
    write ( *, '(i3,3x,i3,3x,i3,3x,i3)' ) xval, less, equal, more
  end do

  return
end
subroutine test0452
!
!*******************************************************************************
!
!! TEST0452 tests IVEC_INDEX_INSERT.
!! TEST0452 tests IVEC_INDEX_DELETE_DUPES.
!! TEST0452 tests IVEC_INDEX_DELETE_ALL.
!! TEST0452 tests IVEC_INDEX_DELETE_ONE.
!
  integer, parameter :: max_n = 25
!
  integer i
  integer indx(max_n)
  integer j
  integer n
  integer n2
  integer x(max_n)
  integer xval
!
  n = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST0452'
  write ( *, * ) '  IVEC_INDEX_INSERT inserts values into an'
  write ( *, * ) '    index sorted array of integers.'
  write ( *, * ) '  IVEC_INDEX_DELETE_ALL deletes all copies of a'
  write ( *, * ) '    particular value.'
  write ( *, * ) '  IVEC_INDEX_DELETE_ONE deletes one copies of a'
  write ( *, * ) '    particular value.'
  write ( *, * ) '  IVEC_INDEX_DELETE_DUPES deletes duplicates.'
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values:'
  write ( *, * ) ' '

  xval = 8
  call ivec_index_insert ( n, x, indx, xval )

  xval = 7
  call ivec_index_insert ( n, x, indx, xval )

  do i = 1, 20
    call i_random ( 0, 20, xval )
    write ( *, '(4x,i3)' ) xval
    call ivec_index_insert ( n, x, indx, xval )
  end do

  xval = 7
  call ivec_index_insert ( n, x, indx, xval )

  xval = 8
  call ivec_index_insert ( n, x, indx, xval )

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Call IVEC_INDEX_DELETE_ONE to delete one value equal to 8:'

  xval = 8
  call ivec_index_delete_one ( n, x, indx, xval )

  write ( *, * ) ' '
  write ( *, * ) '  Call IVEC_INDEX_DELETE_ALL to delete all values equal to 7:'

  xval = 7
  call ivec_index_delete_all ( n, x, indx, xval )

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Call IVEC_INDEX_DELETE_DUPES to delete duplicates:'

  call ivec_index_delete_dupes ( n, x, indx, n2 )

  n = n2

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of unique entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,i3)' ) i, indx(i), x(i)
  end do

  return
end
subroutine test0453
!
!*******************************************************************************
!
!! TEST0453 tests IVEC_INDEX_INSERT_UNIQUE.
!! TEST0453 tests IVEC_INDEX_ORDER.
!
  integer, parameter :: max_n = 20
!
  integer i
  integer indx(max_n)
  integer j
  integer n
  integer x(max_n)
  integer xval
!
  n = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST0453'
  write ( *, * ) '  IVEC_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, * ) '    index sorted array.'
  write ( *, * ) '  IVEC_INDEX_ORDER sorts an index sorted array.'
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values:'
  write ( *, * ) ' '
  do i = 1, 20
    call i_random ( 0, 20, xval )
    write ( *, '(4x,i3)' ) xval
    call ivec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of unique entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,i3,9x,i3)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Now call IVEC_INDEX_ORDER to carry out the sorting:'

  call ivec_index_order ( n, x, indx )

  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  X(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3)' ) i, x(i)
  end do
  return
end
subroutine test046
!
!*******************************************************************************
!
!! TEST046 tests IVEC_MERGE_A;
!! TEST046 tests IVEC_SEARCH_BINARY_A;
!! TEST046 tests IVEC_SORT_HEAP_A.
!
  integer, parameter :: na = 10
  integer, parameter :: nb = 10
!
  integer a(na)
  integer b(nb)
  integer c(na+nb)
  integer index
  integer nc
  integer search_val
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST046'
  write ( *, * ) '  For ascending order:'
  write ( *, * ) '  IVEC_MERGE_A merges two sorted integer arrays;'
  write ( *, * ) '  IVEC_SEARCH_BINARY_A searchs an array for a value;'
  write ( *, * ) '  IVEC_SORT_HEAP_A sorts an integer array.'
  write ( *, * ) ' '
 
  call ivec_random ( 0, na, na, a )
 
  search_val = a(1)

  call ivec_sort_heap_a ( na, a )
 
  call ivec_random ( 0, nb, nb, b )
 
  call ivec_sort_heap_a ( nb, b )
 
  call ivec_print ( na, a, '  Input vector A:' )

  call ivec_print ( nb, b, '  Input vector B:' )
 
  call ivec_merge_a ( na, a, nb, b, nc, c )

  call ivec_print ( nc, c, '  Merged vector C:' )
!
!  Now search the sorted array for a given value.
!
  write ( *, * ) ' '
  write ( *, * ) '  Search the array C for the value ', search_val

  call ivec_search_binary_a ( nc, c, search_val, index )

  write ( *, * ) ' '
  write ( *, * ) '  SEARCH RESULT:'
  if ( index > 0 ) then
    write ( *, * ) '    The value occurs in index ', index
  else
    write ( *, * ) '    The value does not occur in the array.'
  end if

  return
end
subroutine test047
!
!*******************************************************************************
!
!! TEST047 tests IVEC_ORDER_TYPE.
!
  integer, parameter :: n = 4
  integer, parameter :: ntest = 6
!
  integer itest
  integer j
  integer order
  integer x(n,ntest)
!
  x(1,1) = 1
  x(2,1) = 3
  x(3,1) = 2
  x(4,1) = 4

  x(1,2) = 2
  x(2,2) = 2
  x(3,2) = 2
  x(4,2) = 2

  x(1,3) = 1
  x(2,3) = 2
  x(3,3) = 2
  x(4,3) = 4

  x(1,4) = 1
  x(2,4) = 2
  x(3,4) = 3
  x(4,4) = 4

  x(1,5) = 4
  x(2,5) = 4
  x(3,5) = 3
  x(4,5) = 1

  x(1,6) = 9
  x(2,6) = 7
  x(3,6) = 3
  x(4,6) = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST047'
  write ( *, * ) '  IVEC_ORDER_TYPE classifies an integer vector as'
  write ( *, * ) '  -1: no order'
  write ( *, * ) '   0: all equal;'
  write ( *, * ) '   1: ascending;'
  write ( *, * ) '   2: strictly ascending;'
  write ( *, * ) '   3: descending;'
  write ( *, * ) '   4: strictly descending.'
  write ( *, * ) ' '

  do itest = 1, ntest

    call ivec_order_type ( n, x(1,itest), order )

    write ( *, * ) ' '
    write ( *, * ) 'The following vector has order type ', order
    write ( *, * ) ' '
    do j = 1, n
      write ( *, '(i8,i8)' ) j, x(j,itest)
    end do

  end do
   
  return
end
subroutine test0475
!
!*******************************************************************************
!
!! TEST0475 tests IVEC_PAIRWISE_PRIME.
!
  integer, parameter :: n = 4
  integer, parameter :: ntest = 6
!
  integer itest
  logical ivec_pairwise_prime
  integer x(n,ntest)
!
  x(1,1) = 1
  x(2,1) = 3
  x(3,1) = 2
  x(4,1) = 4

  x(1,2) = 2
  x(2,2) = 2
  x(3,2) = 2
  x(4,2) = 2

  x(1,3) = 5
  x(2,3) = 7
  x(3,3) = 12
  x(4,3) = 29

  x(1,4) = 1
  x(2,4) = 13
  x(3,4) = 1
  x(4,4) = 11

  x(1,5) = 1
  x(2,5) = 4
  x(3,5) = 9
  x(4,5) = 16

  x(1,6) = 6
  x(2,6) = 35
  x(3,6) = 13
  x(4,6) = 77

  write ( *, * ) ' '
  write ( *, * ) 'TEST0475'
  write ( *, * ) '  IVEC_PAIRWISE_PRIME determines if a vector of'
  write ( *, * ) '  integers is pairwise prime.'
  write ( *, * ) ' '
  write ( *, * ) '             Pairwise'
  write ( *, * ) 'Row Vector     Prime?'
  write ( *, * ) ' '
  do itest = 1, ntest

    write ( *, '(4i3,3x,l1)' ) x(1:4,itest), ivec_pairwise_prime ( n, x(1,itest) )

  end do
   
  return
end
subroutine test048
!
!*******************************************************************************
!
!! TEST048 tests IVEC_PART.
!
  integer, parameter :: n = 5
!
  integer a(n)
  integer nval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST048'
  write ( *, * ) '  IVEC_PART partitions an integer.'

  nval = 17
  write ( *, * ) ' '
  write ( *, * ) '  NVAL = ', nval

  call ivec_part ( n, nval, a )

  call ivec_print ( n, a, '  Partitioned:' )
 
  nval = - 49
  write ( *, * ) ' '
  write ( *, * ) '  NVAL = ', nval

  call ivec_part ( n, nval, a )

  call ivec_print ( n, a, '  Partitioned:' )

  return
end
subroutine test049
!
!*******************************************************************************
!
!! TEST049 tests IVEC_PART_QUICK_A.
!
  integer, parameter :: n = 12
!
  integer a(n)
  integer i
  integer l
  integer r
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST049'
  write ( *, * ) '  IVEC_PART_QUICK_A reorders an integer vector'
  write ( *, * ) '    as part of a quick sort.'

  call ivec_random ( 0, n, n, a )
 
  call ivec_print ( n, a, '  Before rearrangment:' )

  call ivec_part_quick_a ( n, a, l, r )

  write ( *, * ) ' '
  write ( *, * ) '  Rearranged array'
  write ( *, * ) '  Left = ', l
  write ( *, * ) '  Right = ', r
  write ( *, * ) ' '

  do i = 1, n

    write ( *, '(2i8)' ) i, a(i)

    if ( i == l ) then
      write ( *, * ) ' '
    end if

  end do

  return
end
subroutine test050
!
!*******************************************************************************
!
!! TEST050 tests IVEC_REVERSE.
!
  integer, parameter :: n = 10
!
  integer a(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST050'
  write ( *, * ) '  IVEC_REVERSE reverses a list of integers.'
  write ( *, * ) ' '

  call ivec_random ( 0, 3*n, n, a )
 
  call ivec_print ( n, a, '  Original vector:' )

  call ivec_reverse ( n, a )

  call ivec_print ( n, a, '  Reversed:' )
 
  return
end
subroutine test051
!
!*******************************************************************************
!
!! TEST051 tests IVEC_SORT_HEAP_A.
!! TEST051 tests IVEC_SORT_HEAP_D.
!
  integer, parameter :: n = 20
!
  integer a(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST051'
  write ( *, * ) '  For a vector of integers,'
  write ( *, * ) '  IVEC_SORT_HEAP_A ascending sorts,'
  write ( *, * ) '  IVEC_SORT_HEAP_D descending sorts.'
  write ( *, * ) ' '

  call ivec_random ( 0, 3*n, n, a )

  call ivec_print ( n, a, '  Unsorted:' )

  call ivec_sort_heap_a ( n, a )

  call ivec_print ( n, a, '  Ascending sort:' )
 
  call ivec_sort_heap_d ( n, a )

  call ivec_print ( n, a, '  Descending sort:' )
 
  return
end
subroutine test052
!
!*******************************************************************************
!
!! TEST052 tests IVEC_SORT_HEAP_A;
!! TEST052 tests IVEC_UNIQ.
!
  integer, parameter :: n = 20
!
  integer a(n)
  integer nuniq
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST052'
  write ( *, * ) '  IVEC_SORT_HEAP_A sorts an integer array;'
  write ( *, * ) '  IVEC_UNIQ finds unique entries.'
  write ( *, * ) ' '
 
  call ivec_random ( 0, n, n, a )
 
  call ivec_sort_heap_a ( n, a )
 
  call ivec_print ( n, a, '  Input vector:' )
 
  call ivec_uniq ( n, a, nuniq )
 
  call ivec_print ( nuniq, a, '  Unique entries:' )
 
  return
end
subroutine test053
!
!*******************************************************************************
!
!! TEST053 tests IVEC_SORT_HEAP_INDEX_A.
!! TEST053 tests IVEC_SORT_HEAP_INDEX_D.
!
  integer, parameter :: n = 20
!
  integer a(n)
  integer i
  integer indx(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST053'
  write ( *, * ) '  IVEC_SORT_HEAP_INDEX_A creates an ascending'
  write ( *, * ) '    sort index for an integer array.'
  write ( *, * ) '  IVEC_SORT_HEAP_INDEX_D creates a descending'
  write ( *, * ) '    sort index for an integer array.'

  call ivec_random ( 0, 3 * n, n, a )
 
  call ivec_print ( n, a, '  Unsorted array:' )

  call ivec_sort_heap_index_a ( n, a, indx )

  write ( *, * ) ' '
  write ( *, * ) '  After indexed ascending sort:'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), A(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, indx(i), a(i)
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Now use the index array to carry out the'
  write ( *, * ) '  permutation implicitly.'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), A(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, indx(i), a(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Call IVEC_PERMUTE to carry out the permutation'
  write ( *, * ) '  explicitly.'
  write ( *, * ) ' '

  call ivec_permute ( n, a, indx )

  call ivec_print ( n, a, '  I, A(I)' )

  call ivec_sort_heap_index_d ( n, a, indx )

  write ( *, * ) ' '
  write ( *, * ) '  After indexed descending sort:'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), A(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, indx(i), a(i)
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Now use the index array to carry out the'
  write ( *, * ) '  permutation implicitly.'
  write ( *, * ) ' '
  write ( *, * ) '  INDX(I), A(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(3i8)' ) i, indx(i), a(indx(i))
  end do

  return
end
subroutine test054
!
!*******************************************************************************
!
!! TEST054 tests IVEC_SORT_INSERT_A.
!
  integer, parameter :: n = 10
!
  integer a(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST054'
  write ( *, * ) '  IVEC_SORT_INSERT_A sorts an integer array.'
  write ( *, * ) ' '
 
  call ivec_random ( 0, n, n, a )

  call ivec_print ( n, a, '  Unsorted array:' )

  call ivec_sort_insert_a ( n, a )
 
  call ivec_print ( n, a, '  Sorted array:' )
 
  return
end
subroutine test055
!
!*******************************************************************************
!
!! TEST055 tests IVEC_SORT_QUICK_A.
!
  integer, parameter :: n = 20
!
  integer a(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST055'
  write ( *, * ) '  IVEC_SORT_QUICK_A sorts an integer vector'
  write ( *, * ) '    using quick sort.'
  write ( *, * ) ' '

  call ivec_random ( 0, 3*n, n, a )
 
  call ivec_print ( n, a, '  Unsorted array:' )

  call ivec_sort_quick_a ( n, a )

  call ivec_print ( n, a, '  Sorted array:' )
 
  return
end
subroutine test056
!
!*******************************************************************************
!
!! TEST056 tests IVEC_SORT_SHELL_A.
!
  integer, parameter :: n = 20
!
  integer a(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST056'
  write ( *, * ) '  IVEC_SORT_SHELL_A sorts an integer vector'
  write ( *, * ) '    using Shell''s sort.'
  write ( *, * ) ' '

  call ivec_random ( 0, 3*n, n, a )
 
  call ivec_print ( n, a, '  Unsorted array:' )

  call ivec_sort_shell_a ( n, a )

  call ivec_print ( n, a, '  Sorted array:' )
 
  return
end
subroutine test057
!
!*******************************************************************************
!
!! TEST057 tests IVEC_VALUE_INDEX.
!
  integer, parameter :: max_index = 3
  integer, parameter :: n = 25
!
  integer a(n)
  integer n_index
  integer value
  integer value_index(max_index)
!
  value = 3

  write ( *, * ) ' '
  write ( *, * ) 'TEST057'
  write ( *, * ) '  IVEC_VALUE_INDEX indexes entries equal to'
  write ( *, * ) '  a given value.'
  write ( *, * ) ' '
  write ( *, * ) '  The desired value is ', value
  write ( *, * ) '  Maximum number of indices to find is ', max_index
 
  call ivec_random ( 1, 5, n, a )
 
  call ivec_print ( n, a, '  Input vector A:' )
 
  call ivec_value_index ( n, a, value, max_index, n_index, value_index )
 
  call ivec_print ( n_index, value_index, &
    '  Indices of entries equal to given value: ' )

  return
end
subroutine test058
!
!*******************************************************************************
!
!! TEST058 tests IVEC2_SORT_A.
!! TEST058 tests IVEC2_SORT_D.
!! TEST058 tests IVEC2_UNIQ.
!
  integer, parameter :: n = 10
!
  integer i
  integer ihi
  integer ilo
  integer ivec(n)
  integer jvec(n)
  integer nuniq
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST058'
  write ( *, * ) '  For a pair of integer vectors:'
  write ( *, * ) '  IVEC2_SORT_A ascending sorts;'
  write ( *, * ) '  IVEC2_SORT_D descending sorts;'
  write ( *, * ) '  IVEC2_UNIQ counts unique entries.'
  write ( *, * ) ' '
 
  ilo = 1
  ihi = 3

  call ivec_random ( ilo, ihi, n, ivec )
 
  call ivec_random ( ilo, ihi, n, jvec )

  ivec(3) = ivec(1)
  jvec(3) = jvec(1)

  ivec(5) = ivec(2)
  jvec(5) = jvec(2)

  ivec(9) = ivec(1)
  jvec(9) = jvec(1)

  call ivec2_print ( n, ivec, jvec, '  The array:' )

  call ivec2_sort_a ( n, ivec, jvec )

  call ivec2_print ( n, ivec, jvec, '  After ascending sort:' ) 
 
  call ivec2_sort_d ( n, ivec, jvec )

  call ivec2_print ( n, ivec, jvec, '  After descending sort:' ) 

  call ivec2_uniq ( n, ivec, jvec, nuniq )
 
  call ivec2_print ( nuniq, ivec, jvec, '  After UNIQ:' )

  return
end
subroutine test059
!
!*******************************************************************************
!
!! TEST059 tests NORMAL_01_SAMPLE.
!
  integer i
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST059'
  write ( *, * ) '  NORMAL_01_SAMPLE generates normally distributed'
  write ( *, * ) '    random values.'

  write ( *, * ) ' '

  write ( *, * ) ' '
  do i = 1, 20

    call normal_01_sample ( x )
    write ( *, '(g14.6)' ) x

  end do

  return
end
subroutine test0595
!
!*******************************************************************************
!
!! TEST0595 tests POINTS_NEAREST_POINT_BINS_2D.
!! TEST0595 tests POINTS_NEAREST_POINT_NAIVE_2D.
!
  integer, parameter :: nbin = 20
  integer, parameter :: ndim = 2
  integer, parameter :: nset = 1000
  integer, parameter :: ntest = 10
!
  real, parameter, dimension ( ndim ) :: bin_min = (/  0.0E+00,  0.0E+00 /)
  real, parameter, dimension ( ndim ) :: bin_max = (/ 10.0E+00, 10.0E+00 /)
  integer bin_last(nbin,nbin)
  integer bin_next(nset)
  integer bin_start(nbin,nbin)
  integer compares
  real d
  real d_min1
  real d_min2
  integer i
  integer i_min1
  integer i_min2
  real p(ndim)
  real pset(ndim,nset)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0595'
  write ( *, * ) '  Given a point in 2D, we want to find its nearest'
  write ( *, * ) '  neighbor among points in a set.'
  write ( *, * ) ' '
  write ( *, * ) '  POINTS_NEAREST_POINT_NAIVE_2D uses a naive algorithm.'
  write ( *, * ) '  POINTS_NEAREST_POINT_BINS_2D uses bins.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of points in the pointset is ', nset
  write ( *, * ) '  The number of bins used in each direction is ', nbin
  write ( *, * ) ' '
  write ( *, * ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, * ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, * ) ' '
  write ( *, * ) '  Test point X range:     ', bin_min(1), bin_max(1)
  write ( *, * ) '  Test point Y range:     ', bin_min(2), bin_max(2)
!
!  Set the pointset.
!
  call r2vec_random ( bin_min, bin_max, nset, pset )

! call r2vec_print ( nset, pset, '  The pointset:' )
!
!  Implicitly bin the data.
!
  call r2vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )
!
!  Explicitly reorder the data by bins.
!
  call r2vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )
!
!  Within each bin, sort the data.
!
  call r2vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  Seek nearest neighbors.
!
  write ( *, * ) ' '
  write ( *, * ) '    Test point           Neighbor point      Distance'
  write ( *, * ) '--------------------  --------------------  ----------'
  write ( *, * ) ' '

  do i = 1, ntest

    call r2_random ( bin_min, bin_max, p )

    write ( *, * ) ' '

    call points_nearest_point_naive_2d ( nset, pset, p, i_min1, d_min1 )

    write ( *, '(2f10.4,2x,2f10.4,2x,f10.4)' ) p(1:ndim), pset(1:ndim,i_min1), &
      d_min1

    call points_nearest_point_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
      bin_start, bin_last, bin_next, p, i_min2, d_min2, compares )

    write ( *, '(2f10.4,2x,2f10.4,2x,f10.4,2x,i4)' ) p(1:ndim), &
      pset(1:ndim,i_min2), d_min2, compares

  end do

  return
end
subroutine test0596
!
!*******************************************************************************
!
!! TEST0596 tests POINTS_NEAREST_POINTS_BINS_2D.
!! TEST0596 tests POINTS_NEAREST_POINTS_BINS2_2D.
!! TEST0596 tests POINTS_NEAREST_POINTS_BINS3_2D.
!! TEST0596 tests POINTS_NEAREST_POINTS_NAIVE_2D.
!
  integer, parameter :: nbin = 20
  integer, parameter, dimension ( 2 ) :: nbin_vec = (/ 20, 20 /)
  integer, parameter :: ndim = 2
  integer, parameter :: nset = 1000
  integer, parameter :: ntest = 100
!
  real, parameter, dimension ( ndim ) :: bin_min = (/  0.0E+00,  0.0E+00 /)
  real, parameter, dimension ( ndim ) :: bin_max = (/ 10.0E+00, 10.0E+00 /)
  integer bin_last(nbin,nbin)
  integer bin_next(nset)
  integer bin_start(nbin,nbin)
  integer clock_count1
  integer clock_count2
  integer clock_count3
  integer clock_count4
  integer clock_count5
  integer clock_count6
  integer clock_max
  integer clock_rate
  integer compares(ntest)
  real d
  real d_min0(ntest)
  real d_min1(ntest)
  real d_min2(ntest)
  real d_min3(ntest)
  real eps
  integer i
  integer i_min0(ntest)
  integer i_min1(ntest)
  integer i_min2(ntest)
  integer i_min3(ntest)
  integer n_different
  real pset(ndim,nset)
  real ptest(ndim,ntest)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0596'
  write ( *, * ) '  Given a point set in 2D, and a set of test points,'
  write ( *, * ) '  for each testpoint, find the nearest neighbor in'
  write ( *, * ) '  the point set.'
  write ( *, * ) ' '
  write ( *, * ) '  POINTS_NEAREST_POINTS_NAIVE_2D uses a naive algorithm.'
  write ( *, * ) '  POINTS_NEAREST_POINTS_BINS_2D uses bins.'
  write ( *, * ) '  POINTS_NEAREST_POINTS_BINS2_2D uses bins.'
  write ( *, * ) '  POINTS_NEAREST_POINTS_BINS3_2D uses bins.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of points in the pointset is ', nset
  write ( *, * ) '  The number of bins used in each direction is ', nbin
  write ( *, * ) '  The number of points in the test set is ', ntest
  write ( *, * ) ' '
  write ( *, * ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, * ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, * ) ' '
!
!  Set the pointset.
!
  call r2vec_random ( bin_min, bin_max, nset, pset )
!
!  Set the test points.
!
  call r2vec_random ( bin_min, bin_max, ntest, ptest )
!
!  Implicitly bin the data.
!
  call r2vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )
!
!  Explicitly reorder the data by bins.
!
  call r2vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )
!
!  Within each bin, sort the data.
!
  call r2vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  Seek nearest neighbors.
!
  write ( *, * ) ' '
  write ( *, * ) '  Print results for up to first 10 points...'
  write ( *, * ) ' '
  write ( *, * ) '    Test point                  Distance'
  write ( *, * ) '                       Naive     Bins     Bins2     Bins3'
  write ( *, * ) '--------------------  ------------------------------------'
  write ( *, * ) ' '
  write ( *, * ) ' '

  call system_clock ( clock_count1, clock_rate, clock_max )

  call points_nearest_points_naive_2d ( nset, pset, ntest, ptest, i_min0, &
    d_min0 )

  call system_clock ( clock_count2, clock_rate, clock_max )

  call points_nearest_points_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min1, d_min1, compares )

  call system_clock ( clock_count3, clock_rate, clock_max )

  call points_nearest_points_bins2_2d ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min2, d_min2, compares )

  call system_clock ( clock_count4, clock_rate, clock_max )
!
!  We have to rework the data for BINS3, since we allow a different
!  number of bins in each direction.
!
  call r2vec_bin_even3 ( nset, pset, nbin_vec, bin_min, bin_max, bin_start, &
    bin_last, bin_next )

  call r2vec_binned_reorder ( nset, pset, nbin_vec, bin_start, bin_last, bin_next )

  call r2vec_binned_sort_a ( nset, pset, nbin_vec, bin_start, bin_last )

  call system_clock ( clock_count5, clock_rate, clock_max )

  call points_nearest_points_bins3_2d ( nset, pset, nbin_vec, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min3, d_min3, compares )

  call system_clock ( clock_count6, clock_rate, clock_max )
!
!  Print the results.
!
  do i = 1, min ( ntest, 10 ) 
    write ( *, '(2f10.4,4x,4f10.4)' ) ptest(1:ndim,i), d_min0(i), d_min1(i), &
      d_min2(i), d_min3(i)
  end do
!
!  Check the results.
!
  eps = epsilon ( eps )

  n_different = 0
  do i = 1, ntest
    if ( abs ( d_min0(i) - d_min1(i) ) > eps ) then
      n_different = n_different + 1
    end if
  end do
  
  if ( n_different == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Naive and bin1 codes computed the same results.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) '  Naive and bin1 codes disagreed.'
    write ( *, * ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( abs ( d_min0(i) - d_min2(i) ) > eps ) then
      n_different = n_different + 1
    end if
  end do
  
  if ( n_different == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Naive and bin2 codes computed the same results.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) '  Naive and bin2 codes disagreed.'
    write ( *, * ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( abs ( d_min0(i) - d_min3(i) ) > eps ) then
      n_different = n_different + 1
    end if
  end do
  
  if ( n_different == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Naive and bin3 codes computed the same results.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) '  Naive and bin3 codes disagreed.'
    write ( *, * ) '  Number of discrepancies was ', n_different
  end if

  write ( *, * ) ' '
  write ( *, * ) '  Naive code time = ', &
    real ( clock_count2 - clock_count1 ) / real ( clock_rate )

  write ( *, * ) '  Bin code time =   ', &
    real ( clock_count3 - clock_count2 ) / real ( clock_rate )

  write ( *, * ) '  Bin2 code time =   ', &
    real ( clock_count4 - clock_count3 ) / real ( clock_rate )

  write ( *, * ) '  Bin3 code time =   ', &
    real ( clock_count6 - clock_count5 ) / real ( clock_rate )

  return
end
subroutine test05965
!
!*******************************************************************************
!
!! TEST05965 tests POINTS_NEAREST_POINTS_BINS_2D.
!! TEST05965 tests POINTS_NEAREST_POINTS_BINS2_2D.
!! TEST05965 tests POINTS_NEAREST_POINTS_BINS3_2D.
!! TEST05965 tests POINTS_NEAREST_POINTS_NAIVE_2D.
!
  integer, parameter :: nbin = 10
  integer, parameter, dimension ( 2 ) :: nbin_vec = (/ 4, 25 /)
  integer, parameter :: ndim = 2
  integer, parameter :: nset = 1000
  integer, parameter :: ntest = 100
!
  real, parameter, dimension ( ndim ) :: bin_min = (/  0.0E+00,  0.0E+00 /)
  real, parameter, dimension ( ndim ) :: bin_max = (/ 4.0E+00, 25.0E+00 /)
  integer bin_last(nbin,nbin)
  integer bin_next(nset)
  integer bin_start(nbin,nbin)
  integer clock_count1
  integer clock_count2
  integer clock_count3
  integer clock_count4
  integer clock_count5
  integer clock_count6
  integer clock_max
  integer clock_rate
  integer compares(ntest)
  real d
  real d_min0(ntest)
  real d_min1(ntest)
  real d_min2(ntest)
  real d_min3(ntest)
  real eps
  integer i
  integer i_min0(ntest)
  integer i_min1(ntest)
  integer i_min2(ntest)
  integer i_min3(ntest)
  integer n_different
  real pset(ndim,nset)
  real ptest(ndim,ntest)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05965'
  write ( *, * ) '  Given a point set in 2D, and a set of test points,'
  write ( *, * ) '  for each testpoint, find the nearest neighbor in'
  write ( *, * ) '  the point set.'
  write ( *, * ) ' '
  write ( *, * ) '  In this test, the region is RECTANGULAR.'
  write ( *, * ) '  The BINS and BINS2 codes will end up using rectangular bins;'
  write ( *, * ) '  We will set the BINS3 code to use the same number of bins,'
  write ( *, * ) '  but they will be square.  This should mean that BINS3'
  write ( *, * ) '  finds a match faster.'
  write ( *, * ) ' '
  write ( *, * ) '  POINTS_NEAREST_POINTS_NAIVE_2D uses a naive algorithm.'
  write ( *, * ) '  POINTS_NEAREST_POINTS_BINS_2D uses bins.'
  write ( *, * ) '  POINTS_NEAREST_POINTS_BINS2_2D uses bins.'
  write ( *, * ) '  POINTS_NEAREST_POINTS_BINS3_2D uses bins.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of points in the pointset is ', nset
  write ( *, * ) '  The number of bins used in each direction is ', nbin
  write ( *, * ) '  The number of points in the test set is ', ntest
  write ( *, * ) ' '
  write ( *, * ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, * ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, * ) ' '
!
!  Set the pointset.
!
  call r2vec_random ( bin_min, bin_max, nset, pset )
!
!  Set the test points.
!
  call r2vec_random ( bin_min, bin_max, ntest, ptest )
!
!  Implicitly bin the data.
!
  call r2vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )
!
!  Explicitly reorder the data by bins.
!
  call r2vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )
!
!  Within each bin, sort the data.
!
  call r2vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  Seek nearest neighbors.
!
  write ( *, * ) ' '
  write ( *, * ) '  Print results for up to first 10 points...'
  write ( *, * ) ' '
  write ( *, * ) '    Test point                  Distance'
  write ( *, * ) '                       Naive     Bins     Bins2     Bins3'
  write ( *, * ) '--------------------  ------------------------------------'
  write ( *, * ) ' '
  write ( *, * ) ' '

  call system_clock ( clock_count1, clock_rate, clock_max )

  call points_nearest_points_naive_2d ( nset, pset, ntest, ptest, i_min0, &
    d_min0 )

  call system_clock ( clock_count2, clock_rate, clock_max )

  call points_nearest_points_bins_2d ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min1, d_min1, compares )

  call system_clock ( clock_count3, clock_rate, clock_max )

  call points_nearest_points_bins2_2d ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min2, d_min2, compares )

  call system_clock ( clock_count4, clock_rate, clock_max )
!
!  We have to rework the data for BINS3, since we allow a different
!  number of bins in each direction.
!
  call r2vec_bin_even3 ( nset, pset, nbin_vec, bin_min, bin_max, bin_start, &
    bin_last, bin_next )

  call r2vec_binned_reorder2 ( nset, pset, nbin_vec, bin_start, bin_last, bin_next )

  call r2vec_binned_sort_a2 ( nset, pset, nbin_vec, bin_start, bin_last )

  call system_clock ( clock_count5, clock_rate, clock_max )

  call points_nearest_points_bins3_2d ( nset, pset, nbin_vec, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min3, d_min3, compares )

  call system_clock ( clock_count6, clock_rate, clock_max )
!
!  Print the results.
!
  do i = 1, min ( ntest, 10 )
    write ( *, '(2f10.4,4x,4f10.4)' ) ptest(1:ndim,i), d_min0(i), d_min1(i), &
      d_min2(i), d_min3(i)
  end do
!
!  Check the results.
!
  eps = epsilon ( eps )

  n_different = 0
  do i = 1, ntest
    if ( abs ( d_min0(i) - d_min1(i) ) > eps ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Naive and bin1 codes computed the same results.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) '  Naive and bin1 codes disagreed.'
    write ( *, * ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( abs ( d_min0(i) - d_min2(i) ) > eps ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Naive and bin2 codes computed the same results.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) '  Naive and bin2 codes disagreed.'
    write ( *, * ) '  Number of discrepancies was ', n_different
  end if

  n_different = 0
  do i = 1, ntest
    if ( abs ( d_min0(i) - d_min3(i) ) > eps ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Naive and bin3 codes computed the same results.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) '  Naive and bin3 codes disagreed.'
    write ( *, * ) '  Number of discrepancies was ', n_different
  end if

  write ( *, * ) ' '
  write ( *, * ) '  Naive code time = ', &
    real ( clock_count2 - clock_count1 ) / real ( clock_rate )

  write ( *, * ) '  Bin code time =   ', &
    real ( clock_count3 - clock_count2 ) / real ( clock_rate )

  write ( *, * ) '  Bin2 code time =   ', &
    real ( clock_count4 - clock_count3 ) / real ( clock_rate )

  write ( *, * ) '  Bin3 code time =   ', &
    real ( clock_count6 - clock_count5 ) / real ( clock_rate )

  return
end
subroutine test0597
!
!*******************************************************************************
!
!! TEST0597 tests POINTS_NEAREST_POINTS_BINS2_3D.
!! TEST0597 tests POINTS_NEAREST_POINTS_NAIVE_3D.
!
  integer, parameter :: nbin = 32
  integer, parameter :: ndim = 3
  integer, parameter :: nset = 4096
  integer, parameter :: ntest = 1000
!
  real, parameter, dimension ( ndim ) :: bin_min = (/  0.0E+00,  0.0E+00, 0.0E+00 /)
  real, parameter, dimension ( ndim ) :: bin_max = (/ 10.0E+00, 10.0E+00, 10.0E+00 /)
  integer bin_last(nbin,nbin,nbin)
  integer bin_next(nset)
  integer bin_start(nbin,nbin,nbin)
  integer clock_count1
  integer clock_count2
  integer clock_count3
  integer clock_max
  integer clock_rate
  integer compares(ntest)
  real d
  real d_min1(ntest)
  real d_min2(ntest)
  real eps
  integer i
  integer i_min1(ntest)
  integer i_min2(ntest)
  integer n_different
  real pset(ndim,nset)
  real ptest(ndim,ntest)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0597'
  write ( *, * ) '  Given a point set in 3D, and a set of test points,'
  write ( *, * ) '  for each testpoint, find the nearest neighbor in'
  write ( *, * ) '  the point set.'
  write ( *, * ) ' '
  write ( *, * ) '  POINTS_NEAREST_POINTS_NAIVE_3D uses a naive algorithm.'
  write ( *, * ) '  POINTS_NEAREST_POINTS_BINS2_3D uses bins.'
  write ( *, * ) ' '
  write ( *, * ) '  The number of points in the pointset is ', nset
  write ( *, * ) '  The number of bins used in each direction is ', nbin
  write ( *, * ) '  The number of points in the test set is ', ntest
  write ( *, * ) ' '
  write ( *, * ) '  The X coordinate range: ', bin_min(1), bin_max(1)
  write ( *, * ) '  The Y coordinate range: ', bin_min(2), bin_max(2)
  write ( *, * ) '  The Z coordinate range: ', bin_min(3), bin_max(3)
  write ( *, * ) ' '
!
!  Set the pointset.
!
  call r3vec_random ( bin_min, bin_max, nset, pset )
!
!  Set the test points.
!
  call r3vec_random ( bin_min, bin_max, ntest, ptest )
!
!  Implicitly bin the data.
!
  call r3vec_bin_even2 ( nset, pset, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )
!
!  Explicitly reorder the data by bins.
!
  call r3vec_binned_reorder ( nset, pset, nbin, bin_start, bin_last, bin_next )
!
!  Within each bin, sort the data.
!
  call r3vec_binned_sort_a ( nset, pset, nbin, bin_start, bin_last )
!
!  Seek nearest neighbors.
!
  write ( *, * ) ' '
  write ( *, * ) '  Print up to the first 10 points.'
  write ( *, * ) ' '
  write ( *, * ) '    Test point                       Distance        Comparisons'
  write ( *, * ) '                                 Naive     Bins     Naive Bins'
  write ( *, * ) '-----------------------------  --------------------  ----------'
  write ( *, * ) ' '
  write ( *, * ) ' '

  call system_clock ( clock_count1, clock_rate, clock_max )

  call points_nearest_points_naive_3d ( nset, pset, ntest, ptest, i_min1, &
    d_min1 )

  call system_clock ( clock_count2, clock_rate, clock_max )

  call points_nearest_points_bins2_3d ( nset, pset, nbin, bin_min, bin_max, &
    bin_start, bin_last, bin_next, ntest, ptest, i_min2, d_min2, compares )

  call system_clock ( clock_count3, clock_rate, clock_max )

  do i = 1, min ( ntest, 10 )
    write ( *, '(5f10.4,2x,2i6)' ) ptest(1:ndim,i), d_min1(i), d_min2(i), &
      nset, compares(i)
  end do

  eps = epsilon ( eps )
  n_different = 0
  do i = 1, ntest
    if ( abs ( d_min1(i) - d_min2(i) ) > eps ) then
      n_different = n_different + 1
    end if
  end do

  if ( n_different == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  Naive and bin codes computed the same results.'
  else
    write ( *, * ) ' '
    write ( *, * ) 'WARNING:'
    write ( *, * ) '  Naive and bin codes disagreed.'
    write ( *, * ) '  Number of discrepancies was ', n_different
  end if

  write ( *, * ) ' '
  write ( *, * ) '  Naive code time = ', &
    real ( clock_count2 - clock_count1 ) / real ( clock_rate )

  write ( *, * ) '  Bin code time =   ', &
    real ( clock_count3 - clock_count2 ) / real ( clock_rate )

  return
end
subroutine test060
!
!*******************************************************************************
!
!! TEST060 tests R_DIFF.
!
  integer, parameter :: n = 15
!
  integer i
  integer ndig
  real r_diff
  real y(n)
  real x
!
  y(1) = 0.0625E+00
  y(2) = 0.125E+00
  y(3) = 0.25E+00
  y(4) = 0.50E+00
  y(5) = 0.874E+00
  y(6) = 0.876E+00
  y(7) = 0.90E+00
  y(8) = 0.95E+00
  y(9) = 0.99E+00
  y(10) = 1.0E+00
  y(11) = 1.01E+00
  y(12) = 1.05E+00
  y(13) = 1.10E+00
  y(14) = 3.0E+00
  y(15) = 10.0E+00

  ndig = 3
  x = 1.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST060'
  write ( *, * ) '  R_DIFF computes a difference X-Y to a given'
  write ( *, * ) '    number of binary places.'
  write ( *, * ) ' '
  write ( *, * ) '  For this test, we use ',NDIG,' binary places.'
  write ( *, * ) ' '
  write ( *, * ) '       X       Y       X-Y     R_DIFF(X,Y)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(4f10.4)' ) x, y(i), x-y(i), r_diff(x,y(i),ndig)
  end do
 
  return
end
subroutine test061
!
!*******************************************************************************
!
!! TEST061 tests R_DIGIT.
!
  integer, parameter :: maxdig = 15
!
  integer i
  integer digit(-2:maxdig)
  integer idigit
  real pi
  real x
!
  x = pi ()
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST061'
  write ( *, * ) '  R_DIGIT extracts decimal digits.'
  write ( *, * ) ' '
  write ( *, * ) '  Here, we get digits of ', x
  write ( *, * ) ' '
 
  do idigit = -2, maxdig
    call r_digit ( x, idigit, digit(idigit) )
  end do
 
  write ( *, '(18i3)' ) ( i, i = -2, maxdig )
  write ( *, '(18i3)' ) digit(-2:maxdig)
 
  return
end
subroutine test0615
!
!*******************************************************************************
!
!! TEST0615 tests R_INF.
!
  real r_inf
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0615'
  write ( *, * ) '  R_INF returns real infinity;'
  write ( *, * ) ' '
  write ( *, * ) '  TEST CANCELLED!'
  write ( *, * ) '  I GET A FLOATING-DIVIDE-BY-ZERO ERROR.'
  write ( *, * ) ' '
  return
  write ( *, * ) '    R_INF =   ', r_inf ( )
  write ( *, * ) ' '

  return
end
subroutine test062
!
!*******************************************************************************
!
!! TEST062 tests R_NAN.
!
  real r_nan
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST062'
  write ( *, * ) '  R_NAN returns real NaN.'
  write ( *, * ) ' '
  write ( *, * ) '  TEST CANCELLED!'
  write ( *, * ) '  I GET A FLOATING-DIVIDE-BY-ZERO ERROR.'
  write ( *, * ) ' '
  return

  write ( *, * ) '    R_NAN =   ', r_nan ( )
 
  return
end
subroutine test063
!
!*******************************************************************************
!
!! TEST063 tests R_MANT.
!
  integer is
  integer l
  real r
  real x
!
  x = -314.159E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST063'
  write ( *, * ) '  R_MANT decomposes a real value.'
  write ( *, * ) ' '
  write ( *, * ) '  Number to be decomposed:'
  write ( *, * ) x

  call r_mant ( x, is, r, l )

  write ( *, * ) ' '
  write ( *, * ) '  R_MANT: X = ', is, ' * ', r, ' * 2**', l
 
  return
end
subroutine test0635
!
!*******************************************************************************
!
!! TEST0635 tests R_POWER.
!
  integer i
  integer mults
  integer p
  real r
  real rp
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0635'
  write ( *, * ) '  R_POWER computes R**P, economizing on'
  write ( *, * ) '    multiplications.'
  write ( *, * ) ' '
  write ( *, * ) '      R          P       R**P       Mults'
  write ( *, * ) ' '

  do i = -10, 40

    r = 2.0E+00
    p = i
    call r_power ( r, p, rp, mults )
    write ( *, '(g14.6,i5,g14.6,i5)' ) r, p, rp, mults

  end do

  return
end
subroutine test064
!
!*******************************************************************************
!
!! TEST064 tests R_ROUND2.
!
  integer i
  integer nplace
  real pi
  real x
  real xround
!
  x = pi ( )

  write ( *, * ) ' '
  write ( *, * ) 'TEST064'
  write ( *, * ) '  R_ROUND2 rounds a real number to a'
  write ( *, * ) '    specified number of base 2 digits.'
  write ( *, * ) ' '
  write ( *, * ) '  Test effect on PI:'
  write ( *, * ) '  X = ', x
  write ( *, * ) ' '
  write ( *, * ) '  NPLACE  XROUND'
  write ( *, * ) ' '
 
  do i = 0, 20
    nplace = i
    call r_round2 ( nplace, x, xround )
    write ( *, '(i8,g14.6)' ) i, xround
  end do
 
  return
end
subroutine test065
!
!*******************************************************************************
!
!! TEST065 tests R_ROUNDB.
!
  integer i
  integer ibase
  integer nplace
  real pi
  real x
  real xround
!
  ibase = 3
  x = pi ( )

  write ( *, * ) ' '
  write ( *, * ) 'TEST065'
  write ( *, * ) '  R_ROUNDB rounds a real number to a '
  write ( *, * ) '    specified number of base IBASE digits.'
  write ( *, * ) ' '
  write ( *, * ) '  Here, we will use IBASE = ',ibase
  write ( *, * ) ' '
  write ( *, * ) '  Test effect on PI:'
  write ( *, * ) '  X = ', x
  write ( *, * ) ' '
  write ( *, * ) '  NPLACE  XROUND'
  write ( *, * ) ' '
 
  do i = 0, 20
    nplace = i
    call r_roundb ( ibase, nplace, x, xround )
    write ( *, '(i8,g14.6)' ) i, xround
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Try with a negative base:'
  x = 121.0E+00
  ibase = -3
  nplace = 3
  write ( *, * ) ' '
  write ( *, * ) '  Input quantity is X = ', x
  write ( *, * ) '  to be rounded in base ', ibase
 
  do nplace = 1, 5

    call r_roundb ( ibase, nplace, x, xround )
 
    write ( *, * ) ' '
    write ( *, * ) '  Output value to ', nplace, ' places is ', xround

  end do
 
  return
end
subroutine test066
!
!*******************************************************************************
!
!! TEST066 tests R_ROUNDX.
!
  integer i
  integer nplace
  real pi
  real x
  real xround
!
  x = pi ( )

  write ( *, * ) ' '
  write ( *, * ) 'TEST066'
  write ( *, * ) '  R_ROUNDX rounds a real number to a '
  write ( *, * ) '    specified number of decimal digits.'
  write ( *, * ) ' '
  write ( *, * ) '  Test effect on PI:'
  write ( *, * ) '  X = ', x
  write ( *, * ) ' '
  write ( *, * ) '  NPLACE  XROUND'
  write ( *, * ) ' '
 
  do i = 0, 10
    nplace = i
    call r_roundx ( nplace, x, xround )
    write ( *, '(i8,g16.8)' ) i, xround
  end do
 
  return
end
subroutine test067
!
!*******************************************************************************
!
!! TEST067 tests R_SIGN.
!
  integer, parameter :: test_num = 5
!
  integer test_i
  real r_sign
  real x
  real, parameter, dimension ( test_num ) :: x_test = &
    (/ -1.25E+00, -0.25E+00, 0.0E+00, +0.5E+00, +9.0E+00 /)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST067'
  write ( *, * ) '  R_SIGN returns the sign of a number.'
  write ( *, * ) ' '
 
  do test_i = 1, test_num
    x = x_test(test_i)
    write ( *, '(2f8.4)' ) x, r_sign ( x )
  end do
 
  return
end
subroutine test0675
!
!*******************************************************************************
!
!! TEST0675 tests R_TO_BIN_EVEN.
!! TEST0675 tests BIN_TO_R_EVEN.
!
  real, parameter :: a = 10.0
  real, parameter :: b = 20.0
  integer bin
  real c
  real cmax
  real cmin
  integer i
  integer, parameter :: nbin = 7
  real, parameter :: rmax = 23.0
  real, parameter :: rmin = 8.0
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0675'
  write ( *, * ) '  R_TO_BIN_EVEN puts a number into a bin.'
  write ( *, * ) '  BIN_TO_R_EVEN returns the bin limits.'
  write ( *, * ) '  The bins are equally spaced between A and B,'
  write ( *, * ) '  with two extra bins, for things less than A,'
  write ( *, * ) '  or greater than B.'
  write ( *, * ) ' '
  write ( *, * ) '  A = ', a
  write ( *, * ) '  B = ', b
  write ( *, * ) '  Total number of bins = ', nbin
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values C and put them in bins.'
  write ( *, * ) ' '
  write ( *, * ) '       C      Bin   Bin_Min  Bin_Max'
  write ( *, * ) ' '

  do i = 1, 30
    call r_random ( rmin, rmax, c )
    call r_to_bin_even ( nbin, a, b, c, bin )
    call bin_to_r_even ( nbin, bin, a, b, cmin, cmax )
    write ( *, '(g14.6,i4,2g14.6)' ) c, bin, cmin, cmax
  end do

  return
end
subroutine test0676
!
!*******************************************************************************
!
!! TEST0676 tests R_TO_BIN_EVEN2.
!! TEST0676 tests BIN_TO_R_EVEN2.
!
  real, parameter :: a = 10.0
  real, parameter :: b = 20.0
  integer bin
  real c
  real cmax
  real cmin
  integer i
  integer, parameter :: nbin = 5
  real, parameter :: rmax = 23.0
  real, parameter :: rmin = 8.0
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0676'
  write ( *, * ) '  R_TO_BIN_EVEN2 puts a number into a bin.'
  write ( *, * ) '  BIN_TO_R_EVEN2 returns the bin limits.'
  write ( *, * ) '  The bins are equally spaced between A and B.'
  write ( *, * ) ' '
  write ( *, * ) '  A = ', a
  write ( *, * ) '  B = ', b
  write ( *, * ) '  Total number of bins = ', nbin
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values C and put them in bins.'
  write ( *, * ) ' '
  write ( *, * ) '       C      Bin   Bin_Min  Bin_Max'
  write ( *, * ) ' '

  do i = 1, 30
    call r_random ( rmin, rmax, c )
    call r_to_bin_even2 ( nbin, a, b, c, bin )
    call bin_to_r_even2 ( nbin, bin, a, b, cmin, cmax )
    write ( *, '(g14.6,i4,2g14.6)' ) c, bin, cmin, cmax
  end do

  return
end
subroutine test0677
!
!*******************************************************************************
!
!! TEST0677 tests R2_TO_BIN_EVEN.
!! TEST0677 tests BIN_TO_R2_EVEN.
!
  real, parameter, dimension ( 2 ) :: a = (/  5.0,  0.0 /)
  real, parameter, dimension ( 2 ) :: b = (/ 15.0, 20.0 /)
  integer bin(2)
  real c(2)
  real cmax(2)
  real cmin(2)
  integer i
  integer, parameter :: nbin = 7
  real, parameter, dimension ( 2 ) :: rmin = (/  3.0, -2.0 /)
  real, parameter, dimension ( 2 ) :: rmax = (/ 23.0, 21.0 /)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0677'
  write ( *, * ) '  R2_TO_BIN_EVEN puts an R2 number into a bin.'
  write ( *, * ) '  BIN_TO_R2_EVEN returns the bin limits.'
  write ( *, * ) '  The bins are equally spaced between A and B,'
  write ( *, * ) '  with two extra bins, for things less than A,'
  write ( *, * ) '  or greater than B.'
  write ( *, * ) ' '
  write ( *, * ) '  A(1) = ', a(1)
  write ( *, * ) '  B(1) = ', b(1)
  write ( *, * ) '  A(2) = ', a(2)
  write ( *, * ) '  B(2) = ', b(2)
  write ( *, * ) '  Total number of bins = ', nbin
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values C and put them in bins.'
  write ( *, * ) '  We list the X and Y components on separate lines.'
  write ( *, * ) ' '
  write ( *, * ) '       C      Bin   Bin_Min  Bin_Max'
  write ( *, * ) ' '

  do i = 1, 30
    call r2_random ( rmin, rmax, c )
    call r2_to_bin_even ( nbin, a, b, c, bin )
    call bin_to_r2_even ( nbin, bin, a, b, cmin, cmax )
    write ( *, * ) ' '
    write ( *, '(g14.6,i4,2g14.6)' ) c(1), bin(1), cmin(1), cmax(1)
    write ( *, '(g14.6,i4,2g14.6)' ) c(2), bin(2), cmin(2), cmax(2)
  end do

  return
end
subroutine test068
!
!*******************************************************************************
!
!! TEST068 tests R2_CHEBY.
!
  integer, parameter :: nmax = 10
!
  integer i
  integer n
  real r(nmax)
  real r1
  real r2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST068'
  write ( *, * ) '  R2_CHEBY computes the Chebyshev abscissas'
  write ( *, * ) '    for a given interval [R1,R2].'

  r1 = -1.0E+00
  r2 = +1.0E+00
  n = 5

  call r2_cheby ( n, r, r1, r2 )

  write ( *, * ) ' '
  write ( *, '(a6,g14.6)' ) 'R1:   ', r1
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, r(i)
  end do
  write ( *, '(a6,g14.6)' ) 'R2:   ', r2
 
  r1 =   0.0E+00
  r2 = +10.0E+00
  n = 7

  call r2_cheby ( n, r, r1, r2 )

  write ( *, * ) ' '
  write ( *, '(a6,g14.6)' ) 'R1:   ', r1
  do i = 1, n
    write ( *, '(i6,g14.6)' ) i, r(i)
  end do
  write ( *, '(a6,g14.6)' ) 'R2:   ', r2

  return
end
subroutine test1007
!
!*******************************************************************************
!
!! TEST1007 tests R2VEC_BIN_EVEN.
!! TEST1007 tests R2VEC_BINNED_REORDER.
!! TEST1007 tests R2VEC_BINNED_SORT_A.
!
  integer, parameter :: n = 30
  integer, parameter :: nbin = 4
!
  real a(2,n)
  real, parameter, dimension ( 2 ) :: amin = (/ 8.0E+00, 3.0E+00 /)
  real, parameter, dimension ( 2 ) :: amax = (/ 23.0E+00, 12.0E+00 /)
  real, parameter, dimension ( 2 ) :: bin_min = (/ 10.0E+00, 5.0E+00 /)
  real, parameter, dimension ( 2 ) :: bin_max = (/ 20.0E+00, 10.0E+00 /)
  integer bin_last(nbin,nbin)
  integer bin_next(n)
  integer bin_start(nbin,nbin)
  integer i
  integer i1
  integer i2
  integer j
  integer k
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1007'
  write ( *, * ) '  R2VEC_BIN_EVEN constructs evenly spaced bins and'
  write ( *, * ) '    assigns each element of a real array to a bin.'
  write ( *, * ) '  R2VEC_BINNED_REORDER can reorder the array'
  write ( *, * ) '    to correspond to the bin ordering.'
  write ( *, * ) '  R2VEC_BINNED_SORT_A can sort the individual bins'
  write ( *, * ) '    after the array has been reordered.'
  write ( *, * ) ' '
  write ( *, * ) '  The bins are equally spaced between BIN_MIN and BIN_MAX,'
  write ( *, * ) '  with two extra bins, for things less than BIN_MIN,'
  write ( *, * ) '  or greater than BIN_MAX.'
  write ( *, * ) ' '
  write ( *, * ) '  Component 1 range: ', bin_min(1), bin_max(1)
  write ( *, * ) '  Component 2 range: ', bin_min(2), bin_max(2)
  write ( *, * ) ' '
  write ( *, * ) '  Number of bins per row and column = ', nbin
  write ( *, * ) ' '

  call r2vec_random ( amin, amax, n, a )

  call r2vec_print ( n, a, '  The data vector A to be binned:' )

  call r2vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, bin_last, &
    bin_next )

  call imat_print ( nbin, nbin, nbin, bin_start, '  The BIN_START array:' )

  call imat_print ( nbin, nbin, nbin, bin_start, '  The BIN_LAST array:' )

  call ivec_print ( n, bin_next, '  The BIN_NEXT array:' )

  do i1 = 1, nbin

    do i2 = 1, nbin

      write ( *, * ) ' '
      write ( *, * ) '  Contents of bin number ', i1, i2
      write ( *, * ) ' '

      j = bin_start(i1,i2)
      k = 0

      do while ( j > 0 )
        k = k + 1
        write ( *, '(2i4,2g14.6)' ) k, j, a(1,j), a(2,j)
        j = bin_next(j)
      end do

    end do

  end do
!
!  Now reorder the data to correspond to the bins.
!
  write ( *, * ) ' '
  write ( *, * ) '  Call R2VEC_BINNED_REORDER to reorder the array.'
  write ( *, * ) ' '

  call r2vec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )

  call r2vec_print ( n, a, '  The data vector, sorted by bins:' )

  call imat_print ( nbin, nbin, nbin, bin_start, '  The BIN_START array:' )

  call imat_print ( nbin, nbin, nbin, bin_last, '  The BIN_LAST array:' )

  call ivec_print ( n, bin_next, '  The BIN_NEXT array:' )
!
!  Now sort the bins.
!
  call r2vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

  call r2vec_print ( n, a, '  The data vector, with sorted bins:' )

  return
end
subroutine test0685
!
!*******************************************************************************
!
!! TEST0685 tests RANDOM_SEED.
!! TEST0685 tests RANDOM_NUMBER.
!
  integer i
  integer k
  real mean
  integer, parameter :: n = 1000
  integer seed(1)
  real variance
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0685'
  write ( *, * ) '  RANDOM_SEED is a FORTRAN90 routine which sets or gets'
  write ( *, * ) '    the random number set.'
  write ( *, * ) '  RANDOM_NUMBER returns a uniformly distributed random '
  write ( *, * ) '    value between 0 and 1.'

  k = 1
  call random_seed ( size = k )

  call random_seed ( get = seed )

  write ( *, * ) ' '
  write ( *, * ) '  Starting with seed = ', seed

  call random_number ( x(1:n) )

  call rvec_mean ( n, x, mean )

  call rvec_variance ( n, x, variance )

  write ( *, * ) '  Number of values computed was N = ', n
  write ( *, * ) '  Average value was ', mean
  write ( *, * ) '  Variance was ', variance

  return
end
subroutine test069
!
!*******************************************************************************
!
!! TEST069 tests RAT_FACTOR.
!
  integer, parameter :: maxfactor = 10
!
  integer factor(maxfactor)
  integer i
  integer m
  integer mleft
  integer n
  integer nfactor
  integer nleft
  integer power(maxfactor)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST069'
  write ( *, * ) '  RAT_FACTOR factors a rational value.'

  m = 13 * 7 * 9 * 2
  n = 12
  write ( *, * ) ' '
  write ( *, * ) '  Rational value is ', m, '/', n
  call rat_factor ( m, n, maxfactor, nfactor, factor, power, mleft, nleft )
  write ( *, * ) ' '
  write ( *, * ) '  Prime representation:'
  write ( *, * ) ' '
  write ( *, * ) '  I, FACTOR(I), POWER(I)'
  write ( *, * ) ' '
  if ( mleft /= 1 .or. nleft /= 1 ) then
    write ( *, * ) 0, mleft, ' / ', nleft, ' (UNFACTORED PORTION)'
  end if
  do i = 1, nfactor
    write ( *, '(3i8)' ) i, factor(i), power(i)
  end do
 
  return
end
subroutine test070
!
!*******************************************************************************
!
!! TEST070 tests RCOL_FIND.
!
  integer, parameter :: n = 4
  integer, parameter :: m = 3
  integer, parameter :: lda = m
!
  integer i
  integer icol
  integer j
  integer k
  real rtab(lda,n)
  real rvec(m)

  k = 1
  do i = 1, m
    do j = 1, n

      rtab(i,j) = k

      if ( j == 3 ) then
        rvec(i) = k
      end if

      k = k + 1

    end do
  end do

  call rcol_find ( lda, m, n, rtab, rvec, icol )

  write ( *, * ) ' '
  write ( *, * ) 'TEST070'
  write ( *, * ) '  RCOL_FIND finds a column in a table matching'
  write ( *, * ) '    a given set of data.'
  write ( *, * ) ' '
  write ( *, * ) '  RCOL_FIND returns ICOL = ', icol

  return
end
subroutine test071
!
!*******************************************************************************
!
!! TEST071 tests RCOL_INSERT;
!! TEST071 tests RCOL_SORT_A.
!
  integer, parameter :: maxcol = 10
  integer, parameter :: m = 3
  integer, parameter :: lda = m
!
  real a(lda,maxcol)
  integer i
  integer icol
  integer j
  integer n
  real rvec(m)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST071'
  write ( *, * ) '  RCOL_SORT_A sorts a table of columns.'
  write ( *, * ) '  RCOL_INSERT inserts new columns.'
  write ( *, * ) ' '

  a(1,1) = 2.0E+00
  a(2,1) = 6.0E+00
  a(3,1) = 10.0E+00

  a(1,2) = 4.0E+00
  a(2,2) = 8.0E+00
  a(3,2) = 12.0E+00

  a(1,3) = 1.0E+00
  a(2,3) = 5.0E+00
  a(3,3) = 9.0E+00

  a(1,4) = 3.0E+00
  a(2,4) = 7.0E+00
  a(3,4) = 11.0E+00

  n = 4

  call rmat_print ( lda, m, n, a, '  The unsorted matrix:' )

  call rcol_sort_a ( lda, m, n, a )

  call rmat_print ( lda, m, n, a, '  The sorted matrix:' )

  rvec(1) = 3.0E+00
  rvec(2) = 7.0E+00
  rvec(3) = 11.0E+00

  call rvec_print ( m, rvec, '  New column:' )

  call rcol_insert ( lda, maxcol, m, n, a, rvec, icol )

  if ( icol < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  The data was already in column ', abs ( icol )
  else
    call rmat_print ( lda, m, n, a, '  The updated matrix:' )
  end if

  rvec(1) = 3.0E+00
  rvec(2) = 4.0E+00
  rvec(3) = 18.0E+00

  call rvec_print ( m, rvec, '  New column:' )

  call rcol_insert ( lda, maxcol, m, n, a, rvec, icol )

  if ( icol < 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  The data was already in column ', abs ( icol )
  else
    call rmat_print ( lda, m, n, a, '  The updated matrix:' )
  end if

  return
end
subroutine test072
!
!*******************************************************************************
!
!! TEST072 tests RCOL_MAX;
!! TEST072 tests RCOL_MEAN;
!! TEST072 tests RCOL_MIN;
!! TEST072 tests RCOL_SUM;
!! TEST072 tests RCOL_SWAP;
!! TEST072 tests RCOL_VARIANCE.
!
  integer, parameter :: lda = 5
  integer, parameter :: m = 3
  integer, parameter :: n = 4
!
  real a(lda,n)
  real amax(n)
  real amin(n)
  integer i
  integer icol1
  integer icol2
  integer imax(n)
  integer imin(n)
  integer j
  integer k
  real mean(n)
  real sum(n)
  real variance(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST072'
  write ( *, * ) '  For a real matrix regarded as columns;'
  write ( *, * ) '  RCOL_MAX computes maximums;'
  write ( *, * ) '  RCOL_MEAN computes means;'
  write ( *, * ) '  RCOL_MIN computes minimums;'
  write ( *, * ) '  RCOL_SUM computes sums;'
  write ( *, * ) '  RCOL_SWAP swaps two;'
  write ( *, * ) '  RCOL_VARIANCE computes variances;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k )
    end do
  end do

  call rmat_print ( lda, m, n, a, '  The array:' )

  call rcol_sum ( lda, m, n, a, sum )

  call rcol_max ( lda, m, n, a, imax, amax )

  call rcol_min ( lda, m, n, a, imin, amin )

  call rcol_mean ( lda, m, n, a, mean )

  call rcol_variance ( lda, m, n, a, variance )

  write ( *, * ) ' '
  write ( *, * ) 'Column, maximum, minimum, sum, mean, variance:'
  write ( *, * ) ' '
  do j = 1, n
    write ( *, '(i3,3x,5f10.4)' ) &
      j, amax(j), amin(j), sum(j), mean(j), variance(j)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Swap columns 1 and 3:'

  icol1 = 1
  icol2 = 3
  call rcol_swap ( lda, m, n, a, icol1, icol2 )

  call rmat_print ( lda, m, n, a, '  The updated matrix:' )

  return
end
subroutine test073
!
!*******************************************************************************
!
!! TEST073 tests RCOL_SORTR_A.
!
  integer, parameter :: m = 10
  integer, parameter :: lda = m
  integer, parameter :: n = 3
!
  real a(lda,n)
  integer key
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST073'
  write ( *, * ) '  RCOL_SORTR_A is given an array, and reorders'
  write ( *, * ) '    it so that a particular column is sorted.'
  write ( *, * ) ' '
 
  key = 2
  write ( *, * ) '  Here, the special column is ', key
 
  call rmat_random ( 0.0, 10.0, lda, m, n, a )
 
  call rmat_print ( lda, m, n, a, '  Unsorted array:' )
 
  call rcol_sortr_a ( lda, m, n, a, key )
 
  call rmat_print ( lda, m, n, a, '  Sorted array:' )
 
  return
end
subroutine test1735
!
!*******************************************************************************
!
!! TEST1735 tests RCOL_TO_RVEC.
!
  integer, parameter :: m = 3
  integer, parameter :: n = 4
  integer, parameter :: lda = m
!
  real a(lda,n)
  integer i
  integer j
  real x(m*n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1735'
  write ( *, * ) '  RCOL_TO_RVEC converts an array of columns into a vector.'
  write ( *, * ) ' '
 
  do i = 1, m
    do j = 1, n
      a(i,j) = real ( 10 * i + j )
    end do
  end do

  call rmat_print ( lda, m, n, a, '  The array of columns:' )
 
  call rcol_to_rvec ( lda, m, n, a, x )
 
  call rvec_print ( m*n, x, '  The resulting vector of columns:' )
 
  return
end
subroutine test0775
!
!*******************************************************************************
!
!! TEST0775 tests RMAT_CHOLESKY_FACTOR.
!! TEST0775 tests RMAT_CHOLESKY_SOLVE.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(n)
  real c(lda,n)
  real d(lda,n)
  integer i
  integer ierror
  integer j
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0775'
  write ( *, * ) '  RMAT_CHOLESKY_FACTOR computes the lower'
  write ( *, * ) '  triangular Cholesky factor of a positive'
  write ( *, * ) '  definite symmetric matrix.'
  write ( *, * ) '  RMAT_CHOLESKY_SOLVE solves a linear system'
  write ( *, * ) '  using the Cholesky factorization.'

  do i = 1, n
    do j = 1, n
      if ( i == j ) then
        a(i,j) = 2.0E+00
      else if ( abs ( i - j ) == 1 ) then
        a(i,j) = -1.0E+00
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  call rmat_print ( lda, n, n, a, '  Matrix to be factored:' )

  call rmat_cholesky_factor ( lda, n, a, c, ierror )

  call rmat_print ( lda, n, n, c, '  Cholesky factor C:' )

  call rmat_aat ( lda, lda, n, n, c, d )

  call rmat_print ( lda, n, n, d, '  Product C * C'':' )

  b(1:n-1) = 0.0E+00
  b(n) = real ( n + 1 )

  call rvec_print ( n, b, '  Right hand side:' )

  call rmat_cholesky_solve ( lda, n, c, b, x )

  call rvec_print ( n, x, '  Computed solution:' )

  return
end
subroutine test078
!
!*******************************************************************************
!
!! TEST078 tests RMAT_DET_2D;
!! TEST078 tests RMAT_DET_3D;
!! TEST078 tests RMAT_DET_4D;
!! TEST078 tests RMAT_DET_5D;
!! TEST078 tests RMAT_VAND2.
!
  integer, parameter :: maxn = 5
!
  real a2(2,2)
  real a3(3,3)
  real a4(4,4)
  real a5(5,5)
  real det
  integer i
  integer j
  integer lda
  integer n 
  real rmat_det_2d
  real rmat_det_3d
  real rmat_det_4d
  real rmat_det_5d
  real x(maxn)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST078'
  write ( *, * ) '  Determinants of matrices:'
  write ( *, * ) '  RMAT_DET_2D: 2 by 2 matrix;'
  write ( *, * ) '  RMAT_DET_3D: 3 by 3 matrix;'
  write ( *, * ) '  RMAT_DET_4D: 4 by 4 matrix;'
  write ( *, * ) '  RMAT_DET_5D: 5 by 5 matrix.'

  write ( *, * ) '  RMAT_VAND2 computes the row Vandermonde matrix.'

  x(1) = 1.0E+00
  x(2) = 10.0E+00
  x(3) = 4.0E+00
  x(4) = 2.0E+00
  x(5) = 3.0E+00

  n = 2
  lda = n
  call rmat_vand2 ( lda, n, a2, x )
  det = rmat_det_2d ( a2 )

  call rmat_print ( lda, n, n, a2, '  Matrix:' )

  write ( *, * ) ' '
  write ( *, * ) '  RMAT_DET_2D computes determinant:', det
  det = 1.0E+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, * ) '  Exact determinant is ', det

  n = 3
  lda = n
  call rmat_vand2 ( lda, n, a3, x )
  det = rmat_det_3d ( a3 )

  call rmat_print ( lda, n, n, a3, '  Matrix:' )

  write ( *, * ) ' '
  write ( *, * ) '  RMAT_DET_3D computes determinant:', det
  det = 1.0E+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, * ) '  Exact determinant is ', det

  n = 4
  lda = n
  call rmat_vand2 ( lda, n, a4, x )
  det = rmat_det_4d ( a4 )

  call rmat_print ( lda, n, n, a4, '  Matrix:' )

  write ( *, * ) ' '
  write ( *, * ) '  RMAT_DET_4D computes determinant:', det
  det = 1.0E+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, * ) '  Exact determinant is ', det

  n = 5
  lda = n
  call rmat_vand2 ( lda, n, a5, x )
  det = rmat_det_5d ( a5 )

  call rmat_print ( lda, n, n, a5, '  Matrix:' )

  write ( *, * ) ' '
  write ( *, * ) '  RMAT_DET_5D computes determinant: ', det
  det = 1.0E+00
  do i = 1, n
    do j = 1, i-1
      det = det * ( x(i) - x(j) )
    end do
  end do
  write ( *, * ) '  Exact determinant is ', det

  return
end
subroutine test0745
!
!*******************************************************************************
!
!! TEST0745 tests RMAT_EXPAND_LINEAR.
!
  integer, parameter :: m = 3
  integer, parameter :: n = 2
  integer, parameter :: mb = 10
  integer, parameter :: nb = 5
  integer, parameter :: lda = m
  integer, parameter :: ldb = mb
!
  real a(lda,n)
  real b(ldb,nb)
  integer i
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0745'
  write ( *, * ) '  RMAT_EXPAND_LINEAR fills in a large array by'
  write ( *, * ) '  interpolating data from a small array.'
  write ( *, * ) ' '
 
  do i = 1, m
    do j = 1, n
      a(i,j) = 10.0E+00 * real ( i ) + real ( j )
    end do
  end do

  call rmat_print ( lda, m, n, a, '  The little matrix A:' )
 
  call rmat_expand_linear ( lda, m, n, a, ldb, mb, nb, b )
 
  call rmat_print ( ldb, mb, nb, b, '  Expanded array B:' )
 
  return
end
subroutine test075
!
!*******************************************************************************
!
!! TEST075 tests RMAT_GIVENS_POST.
!
  integer, parameter :: n = 3
  integer, parameter :: lda = n
!
  real a(lda,n)
  real ag(lda,n)
  real g(lda,n)
  integer i
  integer irow
  integer j
  integer jcol
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST075'
  write ( *, * ) '  RMAT_GIVENS_POST computes a Givens ' // &
    'postmultiplier rotation matrices.'
 
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i**(j-1) )
    end do
  end do
 
  call rmat_print ( lda, n, n, a, '  Matrix A:' )
 
  irow = 3
  jcol = 2
 
  write ( *, * ) ' '
  write ( *, * ) '  I, J=', irow, jcol

  call rmat_givens_post ( lda, n, a, g, irow, jcol )
 
  call rmat_print ( lda, n, n, g, '  G' )

  call rmat_mat_mult ( lda, n, a, g, ag )
 
  call rmat_print ( lda, n, n, ag, '  A*G' )

  return
end
subroutine test0756
!
!*******************************************************************************
!
!! TEST0756 tests RMAT_GIVENS_PRE.
!
  integer, parameter :: n = 3
  integer, parameter :: lda = n
!
  real a(lda,n)
  real g(lda,n)
  real ga(lda,n)
  integer i
  integer irow
  integer j
  integer jcol
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0756'
  write ( *, * ) '  RMAT_GIVENS_PRE computes a Givens ' // &
    'premultiplier rotation matrices.'
 
  do i = 1, n
    do j = 1, n
      a(i,j) = real ( i**(j-1) )
    end do
  end do
 
  call rmat_print ( lda, n, n, a, '  Matrix A:' )
 
  irow = 3
  jcol = 2
 
  write ( *, * ) ' '
  write ( *, * ) '  I, J=', irow, jcol
 
  call rmat_givens_pre ( lda, n, a, g, irow, jcol )

  call rmat_print ( lda, n, n, g, '  G' )

  call rmat_mat_mult ( lda, n, g, a, ga )
 
  call rmat_print ( lda, n, n, ga, '  GA' )

  return
end
subroutine test009
!
!*******************************************************************************
!
!! TEST009 tests RVEC_DIF.
!
  integer, parameter :: n = 4
!
  real cof(0:n)
  real fdif
  real :: h = 0.01E+00
  integer i
  real test009_f
  real :: x = 1.0E+00
  real xi
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST009'
  write ( *, * ) '  RVEC_DIF estimates derivatives.'
  write ( *, * ) ' '
  write ( *, * ) '  Estimate the derivative of order N = ', n
  write ( *, * ) '  Using H = ', h
  write ( *, * ) '  at argument X = ', x
!
!  Get the coefficients.
!
  call rvec_dif ( cof, h, n )

  call rvec_print ( n+1, cof, '  The difference coefficients:' )

  fdif = 0.0E+00
  do i = 0, n
    xi = x + real ( 2 * i - n ) * h
    fdif = fdif + cof(i) * test009_f ( xi )
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Estimate is FDIF = ', fdif
 
  return
end
function test009_f ( x )
!
!*******************************************************************************
!
!! TEST009_F evaluates the function used in TEST009.
!
  real test009_f
  real x
!
  test009_f = exp ( x )
 
  return
end
subroutine test0754
!
!*******************************************************************************
!
!! TEST0754 tests RVEC_HOUSE_COLUMN.
!
  integer, parameter :: n = 4
  integer, parameter :: lda = n
!
  real a(lda,n)
  real ahi
  real alo
  real h(lda,n)
  real ha(lda,n)
  integer k
  real v(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0754'
  write ( *, * ) '  RVEC_HOUSE_COLUMN returns the compact form of'
  write ( *, * ) '  a Householder matrix that "packs" a column'
  write ( *, * ) '  of a matrix.'
!
!  Get a random matrix.
!
  alo = 0.0E+00
  ahi = 5.0E+00

  call rmat_random ( alo, ahi, lda, n, n, a )

  call rmat_print ( lda, n, n, a, '  Matrix A:' )

  do k = 1, n-1

    write ( *, * ) ' '
    write ( *, * ) '  Working on column K = ', k

    call rvec_house_column ( n, a(1,k), k, v )

    call rmat_house_form ( lda, n, v, h )

    call rmat_print ( lda, n, n, h, '  Householder matrix H:' )

    call rmat_mat_mult ( lda, n, h, a, ha )

    call rmat_print ( lda, n, n, ha, '  Product H*A:' )
!
!  If we set A := HA, then we can successively convert A to upper
!  triangular form.
!
    a(1:n,1:n) = ha(1:n,1:n)

  end do

  return
end
subroutine test0755
!
!*******************************************************************************
!
!! TEST0755 tests RMAT_HOUSE_FORM.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real h(lda,n)
  real v(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0755'
  write ( *, * ) '  RMAT_HOUSE_FORM forms a Householder'
  write ( *, * ) '    matrix from its compact form.'

  v(1) = 0.0E+00
  v(2) = 0.0E+00
  v(3) = 1.0E+00
  v(4) = 2.0E+00
  v(5) = 3.0E+00

  call rvec_print ( n, v, '  Compact vector form V:' ) 

  call rmat_house_form ( lda, n, v, h )
 
  call rmat_print ( lda, n, n, h, '  Householder matrix H:' )
 
  return
end
subroutine test076
!
!*******************************************************************************
!
!! TEST076 tests RMAT_HOUSE_POST.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real h(lda,n)
  real ha(lda,n)
  integer irow
  integer jcol
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST076'
  write ( *, * ) '  RMAT_HOUSE_POST computes a Householder'
  write ( *, * ) '    postmultiplier;'
 
  call rmat_random ( 0.0E+00, 5.0E+00, lda, n, n, a )
 
  call rmat_print ( lda, n, n, a, '  Matrix A:' )

  irow = 2
  jcol = 3

  write ( *, * ) ' '
  write ( *, * ) '  I, J=', irow, jcol
 
  call rmat_house_post ( lda, n, a, h, irow, jcol )
 
  call rmat_print ( lda, n, n, h, '  Householder matrix H:' )
 
  call rmat_mat_mult ( lda, n, a, h, ha )
 
  call rmat_print ( lda, n, n, ha, '  Product A*H:' )
 
  return
end
subroutine test0761
!
!*******************************************************************************
!
!! TEST0761 tests RMAT_HOUSE_PRE.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real h(lda,n)
  real ha(lda,n)
  integer irow
  integer jcol
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0761'
  write ( *, * ) '  RMAT_HOUSE_PRE computes a Householder'
  write ( *, * ) '    premultiplier;'
 
  call rmat_random ( 0.0E+00, 5.0E+00, lda, n, n, a )
 
  call rmat_print ( lda, n, n, a, '  Matrix A:' )

  irow = 2
  jcol = 3

  write ( *, * ) ' '
  write ( *, * ) '  I, J=', irow, jcol
 
  call rmat_house_pre ( lda, n, a, h, irow, jcol )
 
  call rmat_print ( lda, n, n, h, '  Householder matrix H:' )
 
  call rmat_mat_mult ( lda, n, h, a, ha )
 
  call rmat_print ( lda, n, n, ha, '  Product H*A:' )
 
  return
end
subroutine test082
!
!*******************************************************************************
!
!! TEST082 tests RMAT_INVERSE_2D.
!
  integer, parameter :: n = 2
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  real det
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST082'
  write ( *, * ) '  RMAT_INVERSE_2D inverts a 2 by 2 matrix.'
  write ( *, * ) ' '
!
!  Set the matrix to be inverted.
!
  a(1,1) = 1.0E+00
  a(1,2) = 2.0E+00
 
  a(2,1) = 3.0E+00
  a(2,2) = 4.0E+00

  call rmat_print ( lda, n, n, a, '  Matrix to invert:' )
!
!  Compute the inverse matrix.
!
  call rmat_inverse_2d ( a, b, det )
 
  if ( det == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) '  The input matrix was singular, no inverse'
    write ( *, * ) '  could be computed.'
    write ( *, * ) ' '
    return
  end if

  call rmat_print ( lda, n, n, b, '  Inverse matrix:' )

  call rmat_mat_mult ( lda, n, a, b, c )

  call rmat_print ( lda, n, n, c, '  Product:' )

  return
end
subroutine test083
!
!*******************************************************************************
!
!! TEST083 tests RMAT_INVERSE_3D.
!
  integer, parameter :: n = 3
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  real det
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST083'
  write ( *, * ) '  RMAT_INVERSE_3D inverts a 3 by 3 matrix.'
  write ( *, * ) ' '
!
!  Set the matrix to be inverted.
!
  a(1,1) = 1.0E+00
  a(1,2) = 2.0E+00
  a(1,3) = 3.0E+00
 
  a(2,1) = 4.0E+00
  a(2,2) = 5.0E+00
  a(2,3) = 6.0E+00
 
  a(3,1) = 7.0E+00
  a(3,2) = 8.0E+00
  a(3,3) = 0.0E+00

  call rmat_print ( lda, n, n, a, '  Matrix to be inverted:' )
!
!  Compute the inverse matrix.
!
  call rmat_inverse_3d ( a, b, det )
 
  if ( det == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) '  The input matrix was singular, no inverse'
    write ( *, * ) '  could be computed.'
    write ( *, * ) ' '
    return
  end if

  call rmat_print ( lda, n, n, b, '  Inverse matrix:' )

  call rmat_mat_mult ( lda, n, a, b, c )

  call rmat_print ( lda, n, n, c, '  Product:' )
 
  return
end
subroutine test084
!
!*******************************************************************************
!
!! TEST084 tests RMAT_INVERSE_4D.
!
  integer, parameter :: n = 4
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
  write ( *, * ) 'TEST084'
  write ( *, * ) '  RMAT_INVERSE_4D inverts a 4 x 4 matrix.'

  do i = 1, n
    do j = 1, n

      if ( j >= i ) then
        a(i,j) = real ( n + 1 - j )
      else if ( j == i - 1 ) then
        a(i,j) = n - j
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do

  call rmat_print ( lda, n, n, a, '  Matrix to be inverted:' )

  call rmat_inverse_4d ( a, b, det )

  write ( *, * ) ' '
  write ( *, * ) '  Determinant is ', det

  call rmat_print ( lda, n, n, b, '  Inverse:' )

  call rmat_mat_mult ( lda, n, a, b, c )

  call rmat_print ( lda, n, n, c, '  Product:' )

  return
end
subroutine test0734
!
!*******************************************************************************
!
!! TEST0734 tests RMAT_L_INVERSE.
!
  integer, parameter :: n = 4
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  integer i
  integer j
!
  data ( ( a(i,j), j = 1, n ), i = 1, n ) / &
    1.0E+00, 0.0E+00, 0.0E+00,  0.0E+00, &
    2.0E+00, 3.0E+00, 0.0E+00,  0.0E+00, &
    4.0E+00, 5.0E+00, 6.0E+00,  0.0E+00, &
    7.0E+00, 8.0E+00, 9.0E+00, 10.0E+00 /
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0734'
  write ( *, * ) '  RMAT_L_INVERSE inverts a lower triangular matrix.'

  call rmat_print ( lda, n, n, a, '  Matrix to be inverted:' )
 
  call rmat_l_inverse ( lda, n, a, b )
 
  call rmat_print ( lda, n, n, b, '  Inverse matrix:' )
 
  call rmat_mat_mult ( lda, n, a, b, c )

  call rmat_print ( lda, n, n, c, '  Product:' )

  return
end
subroutine test0735
!
!*******************************************************************************
!
!! TEST0735 tests RMAT_L1_INVERSE.
!
  integer, parameter :: n = 6
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  integer i
  integer j
!
  data ( ( a(i,j), j = 1, n ), i = 1, n ) / &
     1.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, &
     2.0E+00, 1.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, &
     0.0E+00, 0.0E+00, 1.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, &
     5.0E+00, 0.0E+00, 3.0E+00, 1.0E+00, 0.0E+00, 0.0E+00, &
     0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, 1.0E+00, 0.0E+00, &
    75.0E+00, 0.0E+00, 0.0E+00, 6.0E+00, 4.0E+00, 1.0E+00 /
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0735'
  write ( *, * ) '  RMAT_L1_INVERSE inverts a unit lower triangular matrix.'

  call rmat_print ( lda, n, n, a, '  Matrix to be inverted:' )
 
  call rmat_l1_inverse ( lda, n, a, b )
 
  call rmat_print ( lda, n, n, b, '  Inverse matrix:' )
 
  call rmat_mat_mult ( lda, n, a, b, c )

  call rmat_print ( lda, n, n, c, '  Product:' )

  return
end
subroutine test079
!
!*******************************************************************************
!
!! TEST079 tests RMAT_LU.
!
  integer, parameter :: lda = 5
  integer, parameter :: m = 5
  integer, parameter :: n = 5
!
  real a(lda,n)
  integer i
  integer j
  integer k
  integer k2
  real l(lda,m)
  real p(lda,m)
  real prod(lda,n)
  real u(lda,n)
  real, dimension ( n ) :: x = (/ 1.0E+00, 10.0E+00, 4.0E+00, 2.0E+00, &
    3.0E+00 /)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST079'
  write ( *, * ) '  RMAT_LU computes the LU factors of a matrix.'

  call rmat_vand2 ( lda, n, a, x )

  call rmat_print ( lda, m, n, a, '  Matrix to be factored:' )

  call rmat_lu ( lda, m, n, a, l, p, u )

  call rmat_print ( lda, m, m, p, '  P factor:' )

  call rmat_print ( lda, m, m, l, '  L factor:' )
 
  call rmat_print ( lda, m, n, u, '  U factor:' )
 
  do i = 1, m
    do j = 1, n
      prod(i,j) = 0.0E+00
      do k = 1, m
        do k2 = 1, m
          prod(i,j) = prod(i,j) + p(k,i) * l(k,k2) * u(k2,j)
        end do
      end do
    end do
  end do

  call rmat_print ( lda, m, n, prod, '  P*L*U:' )

  return
end
subroutine test074
!
!*******************************************************************************
!
!! TEST074 tests RMAT_IMAX.
!! TEST074 tests RMAT_IMIN.
!! TEST074 tests RMAT_MAX.
!! TEST074 tests RMAT_MAXCOL_MINROW.
!! TEST074 tests RMAT_MAXROW_MINCOL.
!! TEST074 tests RMAT_MIN.
!! TEST074 tests RMAT_MINCOL_MAXROW.
!! TEST074 tests RMAT_MINROW_MAXCOL.
!
  integer, parameter :: m = 5
  integer, parameter :: lda = m
  integer, parameter :: n = 3
!
  real a(lda,n)
  integer i
  integer j
  real rmat_max
  real rmat_maxcol_minrow
  real rmat_maxrow_mincol
  real rmat_min
  real rmat_mincol_maxrow
  real rmat_minrow_maxcol
  real temp1
  real temp2
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST074'
  write ( *, * ) '  RMAT_IMAX locates the maximum;'
  write ( *, * ) '  RMAT_IMIN locates the minimum;'
  write ( *, * ) '  RMAT_MAX computes the maximum;'
  write ( *, * ) '  RMAT_MAXCOL_MINROW computes the maximum over'
  write ( *, * ) '    columns of the mininum over rows;'
  write ( *, * ) '  RMAT_MAXROW_MINCOL computes the maximum over'
  write ( *, * ) '    rows of the mininum over columns;'
  write ( *, * ) '  RMAT_MIN computes the minimum;'
  write ( *, * ) '  RMAT_MINCOL_MAXROW computes the minimum over'
  write ( *, * ) '    columns of the maxinum over rows;'
  write ( *, * ) '  RMAT_MINROW_MAXCOL computes the minimum over'
  write ( *, * ) '    rows of the maxinum over columns;'
  write ( *, * ) ' '
 
  call rmat_random ( 0.0E+00, 10.0E+00, lda, m, n, a )
 
  call rmat_print ( lda, m, n, a, '  Random array:' )
 
  temp1 = rmat_max ( lda, m, n, a )
  temp2 = rmat_min ( lda, m, n, a )

  write ( *, * ) ' '
  write ( *, * ) '  In each pair of numbers, the smaller should'
  write ( *, * ) '  occur first:'
  write ( *, * ) ' '
  write ( *, * ) '  Minimum, Maximum =             ', temp1, temp2
  call rmat_imax ( lda, m, n, a, i, j )
  write ( *, * ) '  Maximum I,J indices            ', i, j
  call rmat_imin ( lda, m, n, a, i, j )
  write ( *, * ) '  Minimum I,J indices            ', i, j

  temp1 = rmat_maxcol_minrow ( lda, m, n, a )
  temp2 = rmat_minrow_maxcol ( lda, m, n, a )

  write ( *, * ) '  MAXCOL_MINROW, MINROW_MAXCOL = ', temp1, temp2

  temp1 = rmat_maxrow_mincol ( lda, m, n, a )
  temp2 = rmat_mincol_maxrow ( lda, m, n, a )

  write ( *, * ) '  MAXROW_MINCOL, MINCOL_MAXROW = ', temp1, temp2

  return
end
subroutine test0799
!
!*******************************************************************************
!
!! TEST0799 tests RMAT_ORTH_RANDOM.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real ata(lda,n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0799'
  write ( *, * ) '  RMAT_ORTH_RANDOM computes a random orthogonal matrix.'

  call rmat_orth_random ( lda, n, a )

  call rmat_print ( lda, n, n, a, '  Random orthogonal matrix A' )

  call rmat_ata ( lda, lda, n, n, a, ata )

  call rmat_print ( lda, n, n, ata, '  AT*A' )

  return
end
subroutine test07995
!
!*******************************************************************************
!
!! TEST07995 tests RMAT_POWER_METHOD.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  real av(n)
  integer i
  integer j
  real r
  real v(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST07995'
  write ( *, * ) '  RMAT_POWER_METHOD applies the power method'
  write ( *, * ) '  to a matrix.'

  do i = 1, n
    do j = 1, n
      if ( j == i - 1 .or. j == i + 1 ) then
        a(i,j) = -1.0E+00
      else if ( j == i ) then
        a(i,j) = 2.0E+00
      else
        a(i,j) = 0.0E+00
      end if
    end do
  end do

  call rmat_power_method ( lda, n, a, r, v )

  write ( *, * ) ' '
  write ( *, * ) '  Estimated eigenvalue = ', r

  call rvec_print ( n, v, '  Estimated eigenvector V:' )

  av(1:n) = matmul ( a(1:n,1:n), v(1:n) )

  call rvec_print ( n, av, '  Value of A*V:' )

  return
end
subroutine test080
!
!*******************************************************************************
!
!! TEST080 tests RMAT_SOLVE.
!
  integer, parameter :: N = 3
  integer, parameter :: NRHS = 2
!
  real a(N,N+NRHS)
  integer i
  integer info
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST080'
  write ( *, * ) '  RMAT_SOLVE solves a 3 by 3 system.'
  write ( *, * ) ' '
!
!  Set the matrix to be inverted.
!
  a(1,1) = 1.0E+00
  a(1,2) = 2.0E+00
  a(1,3) = 3.0E+00
  a(1,4) = 14.0E+00
  a(1,5) = 7.0E+00
 
  a(2,1) = 4.0E+00
  a(2,2) = 5.0E+00
  a(2,3) = 6.0E+00
  a(2,4) = 32.0E+00
  a(2,5) = 16.0E+00
 
  a(3,1) = 7.0E+00
  a(3,2) = 8.0E+00
  a(3,3) = 0.0E+00
  a(3,4) = 23.0E+00
  a(3,5) = 7.0E+00
!
!  Print out the matrix to be inverted.
!
  call rmat_print ( N, N, N+NRHS, a, '  The linear system:' )
!
!  Solve the systems.
!
  call rmat_solve ( a, N, NRHS, info )
 
  if ( info /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) '  The input matrix was singular.'
    write ( *, * ) '  The solutions could not be computed.'
    write ( *, * ) ' '
    return
  else
    write ( *, * ) ' '
    write ( *, * ) '  The computed solutions:'
    write ( *, * ) ' '
    do i = 1, N
      write ( *, '(2g14.6)' ) ( a(i,j), j = N+1, N+NRHS )
    end do
  end if
 
  return
end
subroutine test077
!
!*******************************************************************************
!
!! TEST077 tests RMAT_SYMM_JACOBI;
!! TEST077 tests RMAT_SYMM_RANDOM.
!
  integer, parameter :: n = 5
  integer, parameter :: lda = n
!
  real a(lda,n)
  integer i
  real q(lda,n)
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST077'
  write ( *, * ) '  For a real symmetric matrix:'
  write ( *, * ) '  RMAT_SYMM_JACOBI diagonalizes;'
  write ( *, * ) '  RMAT_SYMM_RANDOM randomizes.'
  write ( *, * ) ' '

  call rvec_identity ( n, x )

  call rmat_symm_random ( lda, n, a, q, x )

  call rmat_print ( lda, n, n, a, '  Matrix to diagonalize:' )

  call rmat_symm_jacobi ( lda, n, a )

  write ( *, * ) ' '
  write ( *, * ) '  Computed Eigenvalues:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(g14.6)' ) a(i,i)
  end do

  return
end
subroutine test081
!
!*******************************************************************************
!
!! TEST081 tests RMAT_TRACE.
!
  integer, parameter :: n = 4
  integer, parameter :: lda = n
!
  real a(lda,n)
  integer i
  integer j
  real rmat_trace
  real trace
!
  do i = 1, n
    do j = 1, n

      if ( j >= i ) then
        a(i,j) = real ( n + 1 - j )
      else if ( j == i - 1 ) then
        a(i,j) = real ( n - j )
      else
        a(i,j) = 0.0E+00
      end if

    end do
  end do
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST081'
  write ( *, * ) '  RMAT_TRACE computes the trace of a matrix'

  call rmat_print ( lda, n, n, a, '  Matrix:' )

  trace = rmat_trace ( lda, n, a )

  write ( *, * ) ' '
  write ( *, * ) '  Trace is ', trace

  return
end
subroutine test0813
!
!*******************************************************************************
!
!! TEST0813 tests RMAT_TRANSPOSE.
!
  integer, parameter :: lda = 6
!
  real a(lda,lda)
  integer i
  integer i_test
  integer j
  integer m
  integer n
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0813'
  write ( *, * ) '  RMAT_TRANSPOSE transposes a matrix'

  do i_test = 1, 3

    if ( i_test == 1 ) then
      m = 4
      n = 2
    else if ( i_test == 2 ) then
      m = 3
      n = 3
    else if ( i_test == 3 ) then
      m = 2
      n = 5
    end if

    do i = 1, m
      do j = 1, n
        a(i,j) = real ( 10 * i + j )
      end do
    end do

    call rmat_print ( lda, m, n, a, '  Matrix:' )

    call rmat_transpose ( lda, m, n, a )

    call i_swap ( m, n )

    call rmat_print ( lda, m, n, a, '  Transposed matrix:' )

  end do

  return
end
subroutine test0814
!
!*******************************************************************************
!
!! TEST0814 tests RMAT_U_INVERSE.
!
  integer, parameter :: n = 4
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  integer i
  integer j
!
  data ( ( a(i,j), j = 1, n ), i = 1, n ) / &
    1.0E+00, 2.0E+00, 4.0E+00,  7.0E+00, &
    0.0E+00, 3.0E+00, 5.0E+00,  8.0E+00, &
    0.0E+00, 0.0E+00, 6.0E+00,  9.0E+00, &
    0.0E+00, 0.0E+00, 0.0E+00, 10.0E+00 /
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0814'
  write ( *, * ) '  RMAT_U_INVERSE inverts an upper triangular matrix.'

  call rmat_print ( lda, n, n, a, '  Input matrix' )
 
  call rmat_u_inverse ( lda, n, a, b )
 
  call rmat_print ( lda, n, n, b, '  Inverse matrix:' )
 
  call rmat_mat_mult ( lda, n, a, b, c )

  call rmat_print ( lda, n, n, c, '  Product:' )

  return
end
subroutine test0815
!
!*******************************************************************************
!
!! TEST815 tests RMAT_U1_INVERSE.
!
  integer, parameter :: n = 6
  integer, parameter :: lda = n
!
  real a(lda,n)
  real b(lda,n)
  real c(lda,n)
  integer i
  integer j
!
  data ( ( a(i,j), j = 1, n ), i = 1, n ) / &
    1.0E+00, 2.0E+00, 0.0E+00, 5.0E+00, 0.0E+00, 75.0E+00, &
    0.0E+00, 1.0E+00, 0.0E+00, 0.0E+00, 0.0E+00,  0.0E+00, &
    0.0E+00, 0.0E+00, 1.0E+00, 3.0E+00, 0.0E+00,  0.0E+00, &
    0.0E+00, 0.0E+00, 0.0E+00, 1.0E+00, 0.0E+00,  6.0E+00, &
    0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, 1.0E+00,  4.0E+00, &
    0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00, 0.0E+00,  1.0E+00 /
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0815'
  write ( *, * ) '  RMAT_U1_INVERSE inverts a unit upper triangular matrix.'

  call rmat_print ( lda, n, n, a, '  Input matrix' )
 
  call rmat_u1_inverse ( lda, n, a, b )
 
  call rmat_print ( lda, n, n, b, '  Inverse matrix:' )

  call rmat_mat_mult ( lda, n, a, b, c )

  call rmat_print ( lda, n, n, c, '  Product:' )

  return
end
subroutine test085
!
!*******************************************************************************
!
!! TEST085 tests ROOTS_TO_RPOLY.
!! TEST085 tests RPOLY_PRINT.
!
  integer, parameter :: n = 4
!
  real c(0:n)
  integer i
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST085'
  write ( *, * ) '  ROOTS_TO_RPOLY computes the coefficients of'
  write ( *, * ) '    a polynomial from its roots.'
  write ( *, * ) '  RPOLY_PRINT prints a polynomial.'

  call rvec_identity ( n, x )

  call rvec_print ( n, x, '  Roots:' )

  call roots_to_rpoly ( n, x, c )

  call rpoly_print ( n, c, '  The polynomial' )

  return
end
subroutine test086
!
!*******************************************************************************
!
!! TEST086 tests RPOLY_LAGRANGE_COEF.
!! TEST086 tests RPOLY_PRINT.
!
  integer, parameter :: npol = 3
!
  integer ipol
  real pcof(0:npol-1)
  real xpol(npol)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST086'
  write ( *, * ) '  RPOLY_LAGRANGE_COEF returns the coefficients'
  write ( *, * ) '    for a Lagrange basis polynomial.'
  write ( *, * ) '  RPOLY_PRINT prints a polynomial.'

  call rvec_identity ( npol, xpol )

  call rvec_print ( npol, xpol, '  Abscissas:' )

  do ipol = 1, npol

    call rpoly_lagrange_coef ( npol, ipol, xpol, pcof )

    call rpoly_print ( npol-1, pcof, '  The Lagrange basis polynomial:' )

  end do

  return
end
subroutine test087
!
!*******************************************************************************
!
!! TEST087 tests RPOLY_LAGRANGE_FACTOR.
!! TEST087 tests RVEC_EVEN.
!
  integer, parameter :: npol = 5
!
  real dwdx
  integer ival
  integer nx
  real wval
  real xhi
  real xlo
  real xpol(npol)
  real xval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST087'
  write ( *, * ) '  RPOLY_LAGRANGE_FACTOR evaluates the Lagrange'
  write ( *, * ) '    factor W(X) at a point.'
  write ( *, * ) ' '
  write ( *, * ) '  For this test, we use ', npol, ' functions.'
!
!  Set the abscissas of the polynomials.
!
  xlo = 0.0E+00
  xhi = real ( npol - 1 )

  call rvec_even ( xlo, xhi, npol, xpol )
 
  call rvec_print ( npol, xpol, '  Abscissas:' )
!
!  Evaluate W(X) and W'(X).
!
  write ( *, * ) ' '
  write ( *, * ) '      X          W(X)          W''(X)'
  write ( *, * ) ' '
 
  nx = 2 * npol - 1
 
  do ival = 1, nx
 
    call rvec_even_select ( xlo, xhi, nx, ival, xval )
 
    call rpoly_lagrange_factor ( npol, xpol, xval, wval, dwdx )
 
    write ( *, '(6g12.4)' ) xval, wval, dwdx
 
  end do
 
  return
end
subroutine test088
!
!*******************************************************************************
!
!! TEST088 tests RPOLY_LAGRANGE_VAL.
!! TEST088 tests RVEC_EVEN.
!
  integer, parameter :: npol = 5
!
  real dpdx(npol)
  integer ipol
  integer ival
  integer nx
  real pval(npol)
  real xhi
  real xlo
  real xpol(npol)
  real xval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST088'
  write ( *, * ) '  RPOLY_LAGRANGE_VAL evaluates a Lagrange'
  write ( *, * ) '    interpolating polynomial at a point.'
  write ( *, * ) ' '
  write ( *, * ) '  For this test, we use ', npol, ' functions.'
!
!  Set the abscissas of the polynomials.
!
  xlo = 0.0E+00
  xhi = real ( npol - 1 )
  call rvec_even ( xlo, xhi, npol, xpol )
 
  call rvec_print ( npol, xpol, '  Abscissas:' )
!
!  Evaluate the polynomials.
!
  write ( *, * ) ' '
  write ( *, * ) '  Here are the values of the functions at '
  write ( *, * ) '  several points:'
  write ( *, * ) ' '
  write ( *, * ) '      X          L1          L2          L3      L4' // &
    '          L5'
  write ( *, * ) ' '
 
  nx = 2 * npol - 1
 
  do ival = 1, nx
 
    call rvec_even_select ( xlo, xhi, nx, ival, xval )
 
    do ipol = 1, npol
      call rpoly_lagrange_val ( npol, ipol, xpol, xval, pval(ipol), dpdx(ipol) )
    end do
 
    write ( *, '(6g12.4)' ) xval, ( pval(ipol), ipol = 1, npol )
 
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  And the derivatives:'
  write ( *, * ) ' '
  write ( *, * ) '      X          L''1         L''2         L''3' // &
    '     L''4         L''5'
  write ( *, * ) ' '
 
  nx = 2 * npol - 1
 
  do ival = 1, nx
 
    call rvec_even_select ( xlo, xhi, nx, ival, xval )
 
    do ipol = 1, npol
      call rpoly_lagrange_val ( npol, ipol, xpol, xval, pval(ipol), dpdx(ipol) )
    end do
 
    write ( *, '(6g12.4)' ) xval, ( dpdx(ipol), ipol = 1, npol )
 
  end do

  return
end
subroutine test089
!
!*******************************************************************************
!
!! TEST089 tests RPOLY_LS_SET;
!! TEST089 tests RPOLY_LS_VAL.
!
  integer, parameter :: npoint = 21
  integer, parameter :: nterms = 4
!
  real b(nterms)
  real c(nterms)
  real d(nterms)
  real f(npoint)
  integer i
  integer k
  real px
  real w(npoint)
  real x(npoint)
!
  w(1:npoint) = 1.0E+00

  do i = 1, npoint
    x(i) = - 1.0E+00 + real ( i - 1 ) / 10.0E+00
    f(i) = real ( int ( exp ( x(i) ) * 100.0E+00 + 0.5E+00 ) ) / 100.0E+00
  end do
 
  call rpoly_ls_set ( b, c, d, f, npoint, nterms, w, x )
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST089'
  write ( *, * ) '  RPOLY_LS_SET sets a least squares polynomial,'
  write ( *, * ) '  RPOLY_LS_VAL evaluates it.'
  write ( *, * ) ' '
  write ( *, * ) '  X, F(X), P(X), Error'
  write ( *, * ) ' '
  do k = 1, nterms
    write ( *, * ) ' '
    write ( *, * ) '  K = ', k
    write ( *, * ) ' '
    do i = 1, npoint
      call rpoly_ls_val ( b, c, d, k, x(i), px )
      write ( *, '(5g14.6)' ) x(i), f(i), px, px - f(i)
    end do
  end do
 
  return
end
subroutine test090
!
!*******************************************************************************
!
!! TEST090 tests RPOLY_LS_SET.
!! TEST090 tests RPOLY_LS_VAL2.
!
  integer, parameter :: npoint = 21
  integer, parameter :: nterms = 4
!
  real b(nterms)
  real c(nterms)
  real d(nterms)
  real f(npoint)
  real fp(npoint)
  integer i
  integer k
  real px
  real pxp
  real w(npoint)
  real x(npoint)
!
  w(1:npoint) = 1.0E+00

  do i = 1, npoint
    x(i) = -1.0E+00 + real ( i-1 ) / 10.0E+00
    f(i) = x(i)**2 - x(i) - 6.0E+00
    fp(i) = 2.0E+00 * x(i) - 1.0E+00
  end do
 
  call rpoly_ls_set ( b, c, d, f, npoint, nterms, w, x )
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST090'
  write ( *, * ) '  RPOLY_LS_SET sets a least squares polynomial,'
  write ( *, * ) '  RPOLY_LS_VAL2 evaluates it.'
  write ( *, * ) ' '
  write ( *, * ) '  X, F(X), P(X), FP(X), PP(X)'
  write ( *, * ) ' '
 
  do k = 1, nterms
    write ( *, * ) ' '
    write ( *, * ) '  K = ', k
    write ( *, * ) ' '
    do i = 1, npoint
      call rpoly_ls_val2 ( b, c, d, k, x(i), px, pxp )
      write ( *, '(5g14.6)' ) x(i), f(i), px, fp(i), pxp
    end do
  end do
 
  return
end
subroutine test091
!
!*******************************************************************************
!
!! TEST091 tests RPOLY_VAL_HORNER.
!
  integer, parameter :: n = 4
!
  integer i
  real, dimension (0:n) :: p = (/ 24.0E+00, -50.0E+00, +35.0E+00, -10.0E+00, &
    1.0E+00 /)
  real pval
  real x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST091'
  write ( *, * ) '  RPOLY_VAL_HORNER evaluates a polynomial at a'
  write ( *, * ) '  point, using Horner''s method.'

  call rpoly_print ( n, p, '  The polynomial:' )

  write ( *, * ) ' '
  write ( *, * ) '  X,  P(X)'
  write ( *, * ) ' '

  do i = 0, 15
    x = real ( i ) / 3.0E+00
    call rpoly_val_horner ( n, p, x, pval )
    write ( *, '(2g14.6)' ) x, pval
  end do
 
  return
end
subroutine test092
!
!*******************************************************************************
!
!! TEST092 tests RPOLY2_EX
!! TEST092 tests RPOLY2_EX2.
!
  real a
  real b
  real c
  integer ierror
  real x1
  real x2
  real x3
  real xmin
  real y1
  real y2
  real y3
  real ymin
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST092'
  write ( *, * ) '  RPOLY2_EX finds the extreme value'
  write ( *, * ) '    of a parabola determined by three points.'
  write ( *, * ) '  RPOLY2_EX2 finds the extreme value'
  write ( *, * ) '    of a parabola determined by three points.'
  a =  2.0E+00
  b = -4.0E+00
  c = 10.0E+00

  x1 = 1.0E+00
  y1 = a * x1**2 + b * x1 + c
  x2 = 2.0E+00
  y2 = a * x2**2 + b * x2 + c
  x3 = 3.0E+00
  y3 = a * x3**2 + b * x3 + c

  write ( *, * ) ' '
  write ( *, * ) '  Parabolic coefficients A, B, C ='
  write ( *, '(3g14.6)' ) a, b, c
  write ( *, * ) ' '
  write ( *, * ) '  X, Y data:'
  write ( *, * ) ' '
  call rvec_print_2d ( x1, y1 )
  call rvec_print_2d ( x2, y2 )
  call rvec_print_2d ( x3, y3 )

  a = 0.0E+00
  b = 0.0E+00
  c = 0.0E+00

  call rpoly2_ex ( x1, y1, x2, y2, x3, y3, xmin, ymin, ierror )

  write ( *, * ) ' '
  write ( *, * ) '  RPOLY2_EX returns XMIN, YMIN = ', xmin, ymin

  call rpoly2_ex2 ( x1, y1, x2, y2, x3, y3, xmin, ymin, a, b, c, ierror )

  write ( *, * ) ' '
  write ( *, * ) '  RPOLY2_EX2 returns XMIN, YMIN = ', xmin, ymin
  write ( *, * ) '  and A, B, C = ', a, b, c

  return
end
subroutine test093
!
!*******************************************************************************
!
!! TEST093 tests RPOLY2_VAL.
!
  integer i
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
  write ( *, * ) ' '
  write ( *, * ) 'TEST093'
  write ( *, * ) '  RPOLY2_VAL evaluates a parabola given'
  write ( *, * ) '    3 data points.'
  write ( *, * ) ' '
  write ( *, * ) '  Our parabola will be 2*x**2 + 3 * x + 1.'
  write ( *, * ) ' '
  write ( *, * ) '  Case 1: 3 distinct data points:'
  write ( *, * ) ' '

  x1 = - 1.0E+00
  x2 = 1.0E+00
  x3 = 3.0E+00

  call test093_f ( x1, y1, yp, ypp )
  call test093_f ( x2, y2, yp, ypp )
  call test093_f ( x3, y3, yp, ypp )

  write ( *, '(5g14.6)' ) x1, y1
  write ( *, '(5g14.6)' ) x2, y2
  write ( *, '(5g14.6)' ) x3, y3

  write ( *, * ) ' '
  write ( *, * ) 'Sampled data:'
  write ( *, * ) ' '
  write ( *, * ) '  X, Y, Y'', Y"'
  write ( *, * ) ' '
  do i = 0, 3
    x = real ( i )
    call rpoly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )
    write ( *, '(4g14.6)' ) x, y, yp, ypp
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Case 2: X1=X2, X3 distinct:'
  write ( *, * ) ' '

  x1 = - 1.0E+00
  x2 = - 1.0E+00
  x3 = 3.0E+00

  call test093_f ( x1, y1, y2, ypp )
  call test093_f ( x3, y3, yp, ypp )

  write ( *, '(2g14.6)' ) x1, y1
  write ( *, '(2g14.6)' ) x2, y2
  write ( *, '(2g14.6)' ) x3, y3

  write ( *, * ) ' '
  write ( *, * ) 'Sampled data:'
  write ( *, * ) ' '
  write ( *, * ) '  X, Y, Y'', Y"'
  write ( *, * ) ' '
  do i = 0, 3
    x = real ( i )
    call rpoly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )
    write ( *, '(4g14.6)' ) x, y, yp, ypp
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Case 3: X1=X2=X3:'
  write ( *, * ) ' '

  x1 = - 1.0E+00
  x2 = - 1.0E+00
  x3 = - 1.0E+00

  call test093_f ( x1, y1, y2, y3 )

  write ( *, '(2g14.6)' ) x1, y1
  write ( *, '(2g14.6)' ) x2, y2
  write ( *, '(2g14.6)' ) x3, y3

  write ( *, * ) ' '
  write ( *, * ) 'Sampled data:'
  write ( *, * ) ' '
  write ( *, * ) '  X, Y, Y'', Y"'
  write ( *, * ) ' '
  do i = 0, 3
    x = real ( i )
    call rpoly2_val ( x1, y1, x2, y2, x3, y3, x, y, yp, ypp )
    write ( *, '(4g14.6)' ) x, y, yp, ypp
  end do

  return
end
subroutine test093_f ( x, y, yp, ypp )
!
!*******************************************************************************
!
!! TEST093_F evaluates a parabola for us.
!
  real x
  real y
  real yp
  real ypp
!
  y = 2.0E+00 * x**2 + 3.0E+00 * x + 1.0E+00
  yp = 4.0E+00 * x + 3.0E+00
  ypp = 4.0E+00

  return
end
subroutine test094
!
!*******************************************************************************
!
!! TEST094 tests RPOLY2_VAL2.
!
  integer, parameter :: ndata = 5
  integer, parameter :: ndim = 2
!
  integer i
  integer left
  real xdata(ndata)
  real xval
  real ydata(ndim,ndata)
  real yval(ndim)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST094'
  write ( *, * ) '  RPOLY2_VAL2 evaluates parabolas through'
  write ( *, * ) '    3 points in a table'
  write ( *, * ) ' '
  write ( *, * ) '  Our data tables will actually be parabolas:'
  write ( *, * ) '    A: 2*x**2 + 3 * x + 1.'
  write ( *, * ) '    B: 4*x**2 - 2 * x + 5.'
  write ( *, * ) ' '

  do i = 1, ndata
    xval = 2.0E+00 * real ( i )
    xdata(i) = xval
    ydata(1,i) = 2.0E+00 * xval**2 + 3.0E+00 * xval + 1.0E+00
    ydata(2,i) = 4.0E+00 * xval**2 - 2.0E+00 * xval + 5.0E+00
    write ( *, '(i6,3g14.6)' ) i, xdata(i), ydata(1,i), ydata(2,i)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Interpolated data:'
  write ( *, * ) ' '
  write ( *, * ) '  LEFT, X, Y1, Y2'
  write ( *, * ) ' '

  do i = 0, 4
    xval = real ( 2 * i + 1 )
    left = max ( min ( i + 1, ndata - 2 ), 1 )
    call rpoly2_val2 ( ndim, ndata, xdata, ydata, left, xval, yval )
    write ( *, '(i8,3g14.6)' ) left, xval, yval(1), yval(2)
  end do

  return
end
subroutine test095
!
!*******************************************************************************
!
!! TEST095 tests RPOLY2_ROOT.
!
  integer, parameter :: ntest = 3
!
  real a(ntest)
  real b(ntest)
  real c(ntest)
  integer i
  complex r1
  complex r2
!
  a(1) = 2.0E+00
  b(1) = -2.0E+00
  c(1) = -24.0E+00
 
  a(2) = 1.0E+00
  b(2) = -20.0E+00
  c(2) = 100.0E+00
 
  a(3) = 1.0E+00
  b(3) = -2.0E+00
  c(3) = 10.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST095'
  write ( *, * ) '  RPOLY2_ROOT finds quadratic equation roots.'
  write ( *, * ) ' '
  write ( *, * ) '         A         B         C     R1         R2'
  write ( *, * ) ' '

  do i = 1, ntest
 
    call rpoly2_root ( a(i), b(i), c(i), r1, r2 )
 
    write ( *, '(3f8.1,4g14.6)' ) a(i), b(i), c(i), r1, r2
 
  end do
 
  return
end
subroutine test096
!
!*******************************************************************************
!
!! TEST096 tests RPOLY3_ROOT.
!
  integer, parameter :: ntest = 4
!
  real a(ntest)
  real b(ntest)
  real c(ntest)
  real d(ntest)
  integer i
  complex r1
  complex r2
  complex r3
!
!  Three distinct real roots, 1, 2, 3.
!
  a(1) = 1.0E+00
  b(1) = -6.0E+00
  c(1) = 11.0E+00
  d(1) = -6.0E+00
!
!  One repeated real root, 1.5, 1.5, 1.5.
!
  a(2) = 8.0E+00
  b(2) = -36.0E+00
  c(2) = 54.0E+00
  d(2) = -27.0E+00
!
!  Two real roots, one repeated, 1, 2, 2.
!
  a(3) = 1.0E+00
  b(3) = -5.0E+00
  c(3) = 8.0E+00
  d(3) = -4.0E+00
!
!  One real root, a complex conjugate pair, 2, 3+2I, 3-2I.
!
  a(4) = 1.0E+00
  b(4) = -8.0E+00
  c(4) = 25.0E+00
  d(4) = -26.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST096'
  write ( *, * ) '  RPOLY3_ROOT finds roots of cubic equations.'
  write ( *, * ) ' '
 
  do i = 1, ntest
 
    write ( *, * ) ' '
    write ( *, * ) '  Polynomial coefficients A, B, C, D:'
    write ( *, * ) ' '
    write ( *, '(4g14.6)' ) a(i), b(i), c(i), d(i)
 
    call rpoly3_root ( a(i), b(i), c(i), d(i), r1, r2, r3 )
 
    write ( *, * ) ' '
    write ( *, * ) '  Roots:'
    write ( *, * ) ' '
    write ( *, '(2g14.6)' ) r1
    write ( *, '(2g14.6)' ) r2
    write ( *, '(2g14.6)' ) r3
 
  end do
 
  return
end
subroutine test097
!
!*******************************************************************************
!
!! TEST097 tests RPOLY4_ROOT.
!
  integer, parameter :: ntest = 7
!
  real a
  real b
  real c
  real d
  real e
  integer i
  complex r1
  complex r2
  complex r3
  complex r4
  real test(5,ntest)
!
!  Four distinct real roots, 1, 2, 3, 4.
!
  test(5,1) = 1.0E+00
  test(4,1) = -10.0E+00
  test(3,1) = +35.0E+00
  test(2,1) = -50.0E+00
  test(1,1) = +24.0E+00
!
!  Three distinct real roots, 1, -2, 3, 3
!
  test(5,2) = 1.0E+00
  test(4,2) = -5.0E+00
  test(3,2) = +1.0E+00
  test(2,2) = +21.0E+00
  test(1,2) = -18.0E+00
!
!  Two distinct real roots, 1, 1, 10, 10.
!
  test(5,3) = 1.0E+00
  test(4,3) = -22.0E+00
  test(3,3) = 141.0E+00
  test(2,3) = -220.0E+00
  test(1,3) = +100.0E+00
!
!  Two distinct real roots, 2, 2, 2, 10
!
  test(5,4) = 1.0E+00
  test(4,4) = -16.0E+00
  test(3,4) = 72.0E+00
  test(2,4) = -128.0E+00
  test(1,4) = 80.0E+00
!
!  One real root, 5, 5, 5, 5
!
  test(5,5) = 1.0E+00
  test(4,5) = -20.0E+00
  test(3,5) = +150.0E+00
  test(2,5) = -500.0E+00
  test(1,5) = +625.0E+00
!
!  Two distinct real roots, one complex conjugate pair.
!
  test(5,6) = 1.0E+00
  test(4,6) = 2.0E+00
  test(3,6) = 1.0E+00
  test(2,6) = 8.0E+00
  test(1,6) = -12.0E+00
!
!  Two distinct complex conjugate pairs.
!
  test(5,7) = 1.0E+00
  test(4,7) = 0.0E+00
  test(3,7) = 13.0E+00
  test(2,7) = 0.0E+00
  test(1,7) = 36.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST097'
  write ( *, * ) '  RPOLY4_ROOT finds roots of quartic equations.'
  write ( *, * ) ' '
 
  do i = 1, ntest
 
    a = test(5,i)
    b = test(4,i)
    c = test(3,i)
    d = test(2,i)
    e = test(1,i)

    write ( *, * ) ' '
    write ( *, * ) '  A, B, C, D, E=', a, b, c, d, e
 
    call rpoly4_root ( a, b, c, d, e, r1, r2, r3, r4 )

    write ( *, * ) ' '
    write ( *, * ) '  Roots:'
    write ( *, * ) ' '
    write ( *, '(5g14.6)' ) r1
    write ( *, '(5g14.6)' ) r2
    write ( *, '(5g14.6)' ) r3
    write ( *, '(5g14.6)' ) r4
 
  end do
 
  return
end
subroutine test098
!
!*******************************************************************************
!
!! TEST098 tests RROW_MAX;
!! TEST098 tests RROW_MEAN;
!! TEST098 tests RROW_MIN;
!! TEST098 tests RROW_SUM;
!! TEST098 tests RROW_SWAP;
!! TEST098 tests RROW_VARIANCE.
!
  integer, parameter :: lda = 5
  integer, parameter :: m = 3
  integer, parameter :: n = 4
!
  real a(lda,n)
  real amax(lda)
  real amin(lda)
  integer i
  integer irow1
  integer irow2
  integer imax(lda)
  integer imin(lda)
  integer j
  integer k
  real mean(m)
  real sum(m)
  real variance(m)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST098'
  write ( *, * ) '  For a real matrix regarded as rows:'
  write ( *, * ) '  RROW_MAX computes maximums;'
  write ( *, * ) '  RROW_MEAN computes means;'
  write ( *, * ) '  RROW_MIN computes minimums;'
  write ( *, * ) '  RROW_SUM computes sums;'
  write ( *, * ) '  RROW_SWAP swaps two;'
  write ( *, * ) '  RROW_VARIANCE computes variances;'

  k = 0
  do i = 1, m
    do j = 1, n
      k = k + 1
      a(i,j) = real ( k )
    end do
  end do

  call rmat_print ( lda, m, n, a, '  The original matrix:' )

  call rrow_sum ( lda, m, n, a, sum )

  call rrow_max ( lda, m, n, a, imax, amax )

  call rrow_min ( lda, m, n, a, imin, amin )

  call rrow_mean ( lda, m, n, a, mean )

  call rrow_variance ( lda, m, n, a, variance )

  write ( *, * ) ' '
  write ( *, * ) 'Maximum, minimum, sum, mean, variance:'
  write ( *, * ) ' '
  do i = 1, m
    write ( *, '(i3,3x,5f10.4)' ) &
      i, amax(i), amin(i), sum(i), mean(i), variance(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) 'Swap rows 1 and 3:'
  write ( *, * ) ' '

  irow1 = 1
  irow2 = 3
  call rrow_swap ( lda, m, n, a, irow1, irow2 )

  call rmat_print ( lda, m, n, a, '  The modified matrix:' )

  return
end
subroutine test0985
!
!*******************************************************************************
!
!! TEST0985 tests RROW_TO_RVEC.
!
  integer, parameter :: m = 3
  integer, parameter :: n = 4
  integer, parameter :: lda = m
!
  real a(lda,n)
  integer i
  integer j
  real x(m*n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST0985'
  write ( *, * ) '  RROW_TO_RVEC converts an array of rows into a vector.'
  write ( *, * ) ' '
 
  do i = 1, m
    do j = 1, n
      a(i,j) = real ( 10 * i + j )
    end do
  end do

  call rmat_print ( lda, m, n, a, '  The array of rows:' )
 
  call rrow_to_rvec ( lda, m, n, a, x )
 
  call rvec_print ( m*n, x, '  The resulting vector of rows:' )
 
  return
end
subroutine test099
!
!*******************************************************************************
!
!! TEST099 tests RVEC_AMAX;
!! TEST099 tests RVEC_AMIN;
!! TEST099 tests RVEC_MAX;
!! TEST099 tests RVEC_MEAN;
!! TEST099 tests RVEC_MEDIAN;
!! TEST099 tests RVEC_MIN.
!! TEST099 tests RVEC_NORM1.
!! TEST099 tests RVEC_NORM2.
!! TEST099 tests RVEC_NORMI.
!! TEST099 tests RVEC_UNIT_SUM.
!
  integer, parameter :: n = 10
!
  real a(n)
  real aval
  integer ival
  real mean
  real median
  real rhi
  real rlo
  real rvec_norm1
  real rvec_norm2
  real rvec_normi
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST099'
  write ( *, * ) '  For a real vector:'
  write ( *, * ) '  RVEC_AMAX:     maximum magnitude entry;'
  write ( *, * ) '  RVEC_AMIN:     minimum magnitude entry.'
  write ( *, * ) '  RVEC_MAX:      maximum entry;'
  write ( *, * ) '  RVEC_MEAN:     mean value;'
  write ( *, * ) '  RVEC_MEDIAN:   median value;'
  write ( *, * ) '  RVEC_MIN:      minimum entry.'
  write ( *, * ) '  RVEC_NORM1:    L1 norm.'
  write ( *, * ) '  RVEC_NORM2:    L2 norm.'
  write ( *, * ) '  RVEC_NORMI:    L-infinity norm.'
  write ( *, * ) '  RVEC_UNIT_SUM: make unit sum;'
  write ( *, * ) ' '
 
  rlo = - real ( n )
  rhi = real ( n )
  call rvec_random ( rlo, rhi, n, a )
 
  call rvec_print ( n, a, '  Input vector:' )

  write ( *, * ) ' '

  call rvec_max ( n, a, aval )
  call rvec_imax ( n, a, ival )
  write ( *, * ) '  Maximum:                 ', aval
  write ( *, * ) '  Maximum index:           ', ival

  call rvec_min ( n, a, aval )
  call rvec_imin ( n, a, ival )
  write ( *, * ) '  Minimum:                 ', aval
  write ( *, * ) '  Minimum index:           ', ival

  call rvec_amax ( n, a, aval )
  call rvec_iamax ( n, a, ival )
  write ( *, * ) '  Maximum absolute:         ', aval
  write ( *, * ) '  Maximum absolute index:   ', ival

  call rvec_amin ( n, a, aval )
  call rvec_iamin ( n, a, ival )
  write ( *, * ) '  Minimum absolute:         ', aval
  write ( *, * ) '  Minimum absolute index:   ', ival

  call rvec_mean ( n, a, mean )
  write ( *, * ) '  Mean:                     ', mean
  call rvec_median ( n, a, median )
  write ( *, * ) '  Median:                   ', median

  call rvec_median ( n, a, median )
  write ( *, * ) '  L1 norm:                 ', rvec_norm1 ( n, a )
  write ( *, * ) '  L2 norm:                 ', rvec_norm2 ( n, a )
  write ( *, * ) '  L-Infinity norm:         ', rvec_normi ( n, a )
 
  call rvec_unit_sum ( n, a )

  call rvec_print ( n, a, '  After calling RVEC_UNIT_SUM:' )

  return
end
subroutine test100
!
!*******************************************************************************
!
!! TEST100 tests RVEC_BIN.
!
  integer, parameter :: n = 25
  integer, parameter :: nbin = 5
!
  integer bin(0:nbin+1)
  real bin_limit(0:nbin)
  real bin_max
  real bin_min
  integer i
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST100'
  write ( *, * ) '  RVEC_BIN computes bins for a real vector.'
 
  call rvec_random ( -2.0E+00, 11.0E+00, n, x )

  call rvec_print ( n, x, '  The vector to be binned:' )

  bin_min =  0.0E+00
  bin_max = 10.0E+00

  write ( *, * ) ' '
  write ( *, * ) '  Number of bins is ', nbin
  write ( *, * ) '  Bin minimum is ', bin_min
  write ( *, * ) '  Bin maximum is ', bin_max

  call rvec_bin ( n, x, nbin, bin_min, bin_max, bin, bin_limit )

  write ( *, * ) ' '
  write ( *, * ) 'Lower Limit    Upper Limit    Count'
  write ( *, * ) ' '

  write ( *, '(2f8.4,i4)' ) bin_min, bin_limit(0), bin(0)
  do i = 1, nbin
    write ( *, '(2f8.4,i4)' ) bin_limit(i-1), bin_limit(i), bin(i)
  end do
  write ( *, '(f8.4,8x,i4)' ) bin_limit(nbin), bin(nbin+1)

  return
end
subroutine test1005
!
!*******************************************************************************
!
!! TEST1005 tests RVEC_BIN_EVEN.
!! TEST1005 tests RVEC_BINNED_REORDER.
!! TEST1005 tests RVEC_BINNED_SORT.
!
  integer, parameter :: n = 30
  integer, parameter :: nbin = 7
!
  real a(n)
  real, parameter :: amax = 23.0E+00
  real, parameter :: amin = 8.0E+00
  real, parameter :: bin_max = 20.0E+00
  real, parameter :: bin_min = 10.0E+00
  integer bin_last(nbin)
  integer bin_next(n)
  integer bin_start(nbin)
  integer i
  integer j
  integer k
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1005'
  write ( *, * ) '  RVEC_BIN_EVEN constructs evenly spaced bins and'
  write ( *, * ) '    assigns each element of a real array to a bin.'
  write ( *, * ) '  RVEC_BINNED_REORDER can reorder the array'
  write ( *, * ) '    to correspond to the bin ordering.'
  write ( *, * ) '  RVEC_BINNED_SORT can sort the array'
  write ( *, * ) '    once it has been reordered.'
  write ( *, * ) ' '
  write ( *, * ) '  The bins are equally spaced between BIN_MIN and BIN_MAX,'
  write ( *, * ) '  with two extra bins, for things less than BIN_MIN,'
  write ( *, * ) '  or greater than BIN_MAX.'
  write ( *, * ) ' '
  write ( *, * ) '  BIN_MIN = ', bin_min
  write ( *, * ) '  BIN_MAX = ', bin_max
  write ( *, * ) '  Total number of bins = ', nbin
  write ( *, * ) ' '

  call rvec_random ( amin, amax, n, a )

  call rvec_print ( n, a, '  The data vector A to be binned:' )

  call rvec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start, &
    bin_last, bin_next )

  call ivec_print ( nbin, bin_start, '  The BIN_START array:' )

  call ivec_print ( nbin, bin_last, '  The BIN_LAST array:' )

  call ivec_print ( n, bin_next, '  The BIN_NEXT array:' )

  do i = 1, nbin

    write ( *, * ) ' '
    write ( *, * ) '  Contents of bin number ', i
    write ( *, * ) ' '

    j = bin_start(i)
    k = 0

    do while ( j > 0 )
      k = k + 1
      write ( *, '(2i4,g14.6)' ) k, j, a(j)
      j = bin_next(j)
    end do

  end do
!
!  Now reorder the data to correspond to the bins.
!
  write ( *, * ) ' '
  write ( *, * ) '  Call RVEC_BINNED_REORDER to reorder the array.'
  write ( *, * ) ' '

  call rvec_binned_reorder ( n, a, nbin, bin_start, bin_last, bin_next )

  call rvec_print ( n, a, '  The data vector A:' )

  call ivec_print ( nbin, bin_start, '  The BIN_START array:' )

  call ivec_print ( nbin, bin_last, '  The BIN_LAST array:' )

  call ivec_print ( n, bin_next, '  The BIN_NEXT array:' )
!
!  Now sort the data, one bin at a time
!
  call rvec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

  call rvec_print ( n, a, '  The sorted data vector A:' )

  return
end
subroutine test101
!
!*******************************************************************************
!
!! TEST101 tests RVEC_BRACKET.
!
  integer, parameter :: n = 10
  integer, parameter :: ntest = 6
!
  integer i
  integer itest
  integer left
  integer right
  real x(n)
  real xtest(ntest)
  real xval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST101'
  write ( *, * ) '  RVEC_BRACKET finds a pair of entries in a'
  write ( *, * ) '    sorted real array which bracket a value.'

  xtest(1) = -10.0E+00
  xtest(2) = 1.0E+00
  xtest(3) = 4.5E+00
  xtest(4) = 5.0E+00
  xtest(5) = 10.0E+00
  xtest(6) = 12.0E+00

  call rvec_identity ( n, x )
  x(6) = x(5)

  call rvec_print ( n, x, '  Sorted array:' )

  write ( *, * ) ' '
  write ( *, * ) '    LEFT             RIGHT'
  write ( *, * ) '  X(LEFT)   XVAL   X(RIGHT)'
  write ( *, * ) ' '

  do itest = 1, ntest

    xval = xtest(itest)

    call rvec_bracket ( n, x, xval, left, right )

    write ( *, '(i14,14x,i14)' ) left, right

    if ( left >= 1 .and. right >= 1 ) then
      write ( *, '(3g14.6)' ) x(left), xval, x(right)
    else if ( left < 1 .and. right >= 1 ) then
      write ( *, '(14x,2g14.6)' )          xval, x(right)
    else if ( left >= 1 .and. right < 1 ) then
      write ( *, '(2g14.6)' ) x(left), xval
    else if ( left < 1 .and. right < 1 ) then
      write ( *, '(14x,g14.6)' )          xval
    end if

  end do

  return
end
subroutine test102
!
!*******************************************************************************
!
!! TEST102 tests RVEC_BRACKET2.
!
  integer, parameter :: n = 10
  integer, parameter :: ntest = 6
!
  integer i
  integer itest
  integer left
  integer right
  integer start
  real x(n)
  real xtest(ntest)
  real xval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST102'
  write ( *, * ) '  RVEC_BRACKET2 finds a pair of entries in a'
  write ( *, * ) '    sorted real array which bracket a value.'

  xtest(1) = -10.0E+00
  xtest(2) = 1.0E+00
  xtest(3) = 4.5E+00
  xtest(4) = 5.0E+00
  xtest(5) = 10.0E+00
  xtest(6) = 12.0E+00

  call rvec_identity ( n, x )
  x(6) = x(5)

  call rvec_print ( n, x, '  Sorted array:' )

  left = 0

  do itest = 1, ntest

    xval = xtest(itest)

    write ( *, * ) ' '
    write ( *, * ) 'Search for XVAL = ', xval

    if ( left > 0 ) then
      start = left
    else
      start = ( n + 1 ) / 2
    end if

    write ( *, * ) 'Start = ', start

    call rvec_bracket2 ( n, x, xval, start, left, right )

    write ( *, * ) 'Left = ', left
    write ( *, * ) 'Right = ', right

    if ( left >= 1 ) then 
      write ( *, * ) 'X(LEFT)=', x(left)
    end if

    if ( right >= 1 ) then
      write ( *, * ) 'X(RIGHT) = ', x(right)
    end if

  end do

  return
end
subroutine test103
!
!*******************************************************************************
!
!! TEST103 tests RVEC_BRACKET3.
!
  integer, parameter :: n = 10
  integer, parameter :: ntest = 6
!
  integer i
  integer itest
  integer left
  real x(n)
  real xtest(ntest)
  real xval
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST103'
  write ( *, * ) '  RVEC_BRACKET3 finds a pair of entries in a'
  write ( *, * ) '    sorted real array which bracket a value.'

  xtest(1) = -10.0E+00
  xtest(2) = 1.0E+00
  xtest(3) = 4.5E+00
  xtest(4) = 5.0E+00
  xtest(5) = 10.0E+00
  xtest(6) = 12.0E+00

  call rvec_identity ( n, x )
  x(6) = x(5)

  call rvec_print ( n, x, '  Sorted array:' )

  left = ( n + 1 ) / 2

  do itest = 1, ntest

    xval = xtest(itest)

    write ( *, * ) ' '
    write ( *, * ) '  Search for XVAL = ', xval

    write ( *, * ) '  Starting guess for interval is = ', left

    call rvec_bracket3 ( n, x, xval, left )

    write ( *, * ) '  Nearest interval:'
    write ( *, * ) '    X{', left,' ]= ', x(left)
    write ( *, * ) '    X[', left+1, ' ]= ', x(left+1)

  end do

  return
end
subroutine test1035
!
!*******************************************************************************
!
!! TEST1035 tests RVEC_CONVOLVE_CIRC
!
  integer, parameter :: n = 4
!
  real x(n)
  real y(n)
  real z(n)
  real z_correct(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1035'
  write ( *, * ) '  RVEC_CONVOLVE_CIRC computes the circular convolution'
  write ( *, * ) '  of two vectors.'

  x(1:4) = (/ 1.0E+00, 2.0E+00, 3.0E+00, 4.0E+00 /)
  y(1:4) = (/ 1.0E+00, 2.0E+00, 4.0E+00, 8.0E+00 /)

  call rvec_print ( n, x, '  The factor X:' )
  call rvec_print ( n, y, '  The factor Y:' )

  call rvec_convolve_circ ( n, x, y, z )

  call rvec_print ( n, z, '  The circular convolution z = xCCy:' )

  z_correct(1:4) = (/ 37.0E+00, 44.0E+00, 43.0E+00, 26.0E+00 /)

  call rvec_print ( n, z_correct, '  Correct answer:' )

  return
end
subroutine test104
!
!*******************************************************************************
!
!! TEST104 tests RVEC_EVEN.
!
  integer, parameter :: n = 10
!
  real x(n)
  real xhi
  real xlo
!
  xlo = 0.0E+00
  xhi = 99.0E+00
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST104'
  write ( *, * ) '  RVEC_EVEN computes N evenly spaced values'
  write ( *, * ) '    between XLO and XHI.'
  write ( *, * ) ' '
  write ( *, * ) '  XLO = ', xlo, ' and XHI = ', xhi
  write ( *, * ) '  while N = ', n
 
  call rvec_even ( xlo, xhi, n, x )
 
  call rvec_print ( n, x, '  Resulting array:' )
 
  return
end
subroutine test105
!
!*******************************************************************************
!
!! TEST105 tests RVEC_EVEN2.
!
  integer, parameter :: nold = 5
  integer, parameter :: maxval = 20
!
  integer i
  integer istar
  integer jstar
  integer nfill(nold-1)
  integer nval
  real xold(nold)
  real xval(maxval)
!
  xold(1) = 0.0E+00
  nfill(1) = 4
  xold(2) = 1.0E+00
  nfill(2) = 3
  xold(3) = 5.0E+00
  nfill(3) = 5
  xold(4) = 2.0E+00
  nfill(4) = 0
  xold(5) = 0.0E+00
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST105'
  write ( *, * ) '  RVEC_EVEN2 interpolates a specified number of '
  write ( *, * ) '    points pairs of values in a vector.'
  write ( *, * ) ' '
  write ( *, * ) '  Input data:'
  write ( *, * ) ' '
  do i = 1, nold
    write ( *, '(g14.6)' ) xold(i)
    if ( i < nold ) then
      write ( *, '(i10)' ) nfill(i)
    end if
  end do
 
  call rvec_even2 ( maxval, nfill, nold, nval, xold, xval )
 
  write ( *, * ) ' '
  write ( *, * ) '  Resulting vector:'
  write ( *, * ) ' '
 
  istar = 1
  jstar = 1
  do i = 1, nval
 
    if ( i == istar ) then
 
      write ( *, '(a1,g14.6)' ) '*', xval(i)
 
      if ( jstar < nold ) then
        istar = istar + nfill(jstar) + 1
        jstar = jstar + 1
      end if
 
    else
 
      write ( *, '(g14.6)' ) xval(i)
 
    end if
 
  end do
 
  return
end
subroutine test106
!
!*******************************************************************************
!
!! TEST106 tests RVEC_EVEN3.
!
  integer, parameter :: nold = 4
  integer, parameter :: nval = 12
!
  real xold(nold)
  real xval(nval)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST106'
  write ( *, * ) '  RVEC_EVEN3 tries to evenly interpolate new data'
  write ( *, * ) '    between old values.'
  write ( *, * ) ' '

  xold(1) = 0.0E+00
  xold(2) = 5.1E+00
  xold(3) = 7.0E+00
  xold(4) = 10.0E+00
 
  call rvec_print ( nold, xold, '  Original vector:' )

  call rvec_even3 ( nold, nval, xold, xval )
 
  call rvec_print ( nval, xval, '  New vector:' )
 
  return
end
subroutine test107
!
!*******************************************************************************
!
!! TEST107 tests RVEC_FRAC.
!
  integer, parameter :: n = 10
!
  real a(n)
  real afrac
  integer k
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST107'
  write ( *, * ) '  RVEC_FRAC: K-th smallest real vector entry;'
  write ( *, * ) ' '

  call rvec_random ( 0.0E+00, 1.0E+00, n, a )

  call rvec_print ( n, a, '  Array to search:' )

  write ( *, * ) ' '
  write ( *, * ) 'Fractile  Value '
  write ( *, * ) ' '

  do k = 1, n, n/2

    call rvec_frac ( n, a, k, afrac )

    write ( *, '(i8,g14.6)' ) k, afrac

  end do

  return
end
subroutine test1071
!
!*******************************************************************************
!
!! TEST1071 tests RVEC_INDEX_SEARCH.
!! TEST1071 tests RVEC_INDEX_INSERT_UNIQUE.
!
  integer, parameter :: max_n = 20
!
  integer equal
  integer i
  integer indx(max_n)
  integer j
  integer less
  integer more
  integer n
  real x(max_n)
  real xval
!
  n = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST1071'
  write ( *, * ) '  RVEC_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, * ) '    index sorted array.'
  write ( *, * ) '  RVEC_INDEX_SEARCH searches for an entry with a given value.'
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values:'
  write ( *, * ) ' '
  do i = 1, 20
    call r_random ( 0.0E+00, 20.0E+00, xval )
    xval = real ( nint ( xval ) )
    write ( *, '(4x,f6.2)' ) xval
    call rvec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Results of search for given XVAL:'
  write ( *, * ) ' '
  write ( *, '(a)' ) 'XVAL  Less Equal More'
  write ( *, * ) ' '

  do i = 0, 20
    xval = real ( i )
    call rvec_index_search ( n, x, indx, xval, less, equal, more )
    write ( *, '(f6.2,3x,i3,3x,i3,3x,i3)' ) xval, less, equal, more
  end do

  return
end
subroutine test1072
!
!*******************************************************************************
!
!! TEST1072 tests RVEC_INDEX_INSERT.
!! TEST1072 tests RVEC_INDEX_DELETE_DUPES.
!! TEST1072 tests RVEC_INDEX_DELETE_ALL.
!! TEST1072 tests RVEC_INDEX_DELETE_ONE.
!
  integer, parameter :: max_n = 25
!
  integer i
  integer indx(max_n)
  integer j
  integer n
  integer n2
  real x(max_n)
  real xval
!
  n = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST1072'
  write ( *, * ) '  RVEC_INDEX_INSERT inserts values into an'
  write ( *, * ) '    index sorted array.'
  write ( *, * ) '  RVEC_INDEX_DELETE_ALL deletes all copies of a'
  write ( *, * ) '    particular value.'
  write ( *, * ) '  RVEC_INDEX_DELETE_ONE deletes one copies of a'
  write ( *, * ) '    particular value.'
  write ( *, * ) '  RVEC_INDEX_DELETE_DUPES deletes duplicates.'
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values:'
  write ( *, * ) ' '

  xval = 8.0E+00
  call rvec_index_insert ( n, x, indx, xval )

  xval = 7.0E+00
  call rvec_index_insert ( n, x, indx, xval )

  do i = 1, 20
    call r_random ( 0.0E+00, 20.0E+00, xval )
    xval = real ( nint ( xval ) )
    write ( *, '(4x,f6.2)' ) xval
    call rvec_index_insert ( n, x, indx, xval )
  end do

  xval = 7.0E+00
  call rvec_index_insert ( n, x, indx, xval )

  xval = 8.0E+00
  call rvec_index_insert ( n, x, indx, xval )

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Call RVEC_INDEX_DELETE_ONE to delete one value equal to 8:'

  xval = 8.0E+00
  call rvec_index_delete_one ( n, x, indx, xval )

  write ( *, * ) ' '
  write ( *, * ) '  Call RVEC_INDEX_DELETE_ALL to delete all values equal to 7:'

  xval = 7.0E+00
  call rvec_index_delete_all ( n, x, indx, xval )

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Call RVEC_INDEX_DELETE_DUPES to delete duplicates:'

  call rvec_index_delete_dupes ( n, x, indx, n2 )

  n = n2

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of unique entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,f6.2)' ) i, indx(i), x(i)
  end do

  return
end
subroutine test1073
!
!*******************************************************************************
!
!! TEST1073 tests RVEC_INDEX_INSERT_UNIQUE.
!! TEST1073 tests RVEC_INDEX_ORDER.
!
  integer, parameter :: max_n = 20
!
  integer i
  integer indx(max_n)
  integer j
  integer n
  real x(max_n)
  real xval
!
  n = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST1073'
  write ( *, * ) '  RVEC_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, * ) '    index sorted array.'
  write ( *, * ) '  RVEC_INDEX_ORDER sorts an index sorted array.'
  write ( *, * ) ' '
  write ( *, * ) '  Generate some random values:'
  write ( *, * ) ' '
  do i = 1, 20
    call r_random ( 0.0E+00, 20.0E+00, xval )
    xval = real ( nint ( xval ) )
    write ( *, '(4x,f6.2)' ) xval
    call rvec_index_insert_unique ( n, x, indx, xval )
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of unique entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,f6.2,9x,f6.2)' ) i, indx(i), x(i), x(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Now call RVEC_INDEX_ORDER to carry out the sorting:'

  call rvec_index_order ( n, x, indx )

  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  X(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,f6.2)' ) i, x(i)
  end do
  return
end
subroutine test108
!
!*******************************************************************************
!
!! TEST108 tests RVEC_MERGE_A;
!! TEST108 tests RVEC_SEARCH_BINARY_A;
!! TEST108 tests RVEC_SORT_HEAP_A.
!
  integer, parameter :: na = 10
  integer, parameter :: nb = 10
!
  real a(na)
  real b(nb)
  real c(na+nb)
  integer index
  integer nc
  real search_val
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST108'
  write ( *, * ) '  For ascending order:'
  write ( *, * ) '  RVEC_MERGE_A merges two sorted real arrays;'
  write ( *, * ) '  RVEC_SEARCH_BINARY_A searches a sorted array;'
  write ( *, * ) '  RVEC_SORT_HEAP_A sorts a real array.'
  write ( *, * ) ' '

  call rvec_random ( 0.0E+00, 1.0E+00, na, a )
  call rvec_random ( 0.0E+00, 1.0E+00, nb, b )

  search_val = a(1)
 
  call rvec_sort_heap_a ( na, a )
 
  call rvec_sort_heap_a ( nb, b )

  call rvec_print ( na, a, '  Sorted vector A:' )

  call rvec_print ( nb, b, '  Sorted vector B:' )
 
  call rvec_merge_a ( na, a, nb, b, nc, c )

  call rvec_print ( nc, c, '  Merged vector C:' )
!
!  Now search the sorted array for a given value.
!
  write ( *, * ) ' '
  write ( *, * ) '  Search the array C for the value ', search_val

  call rvec_search_binary_a ( nc, c, search_val, index )

  write ( *, * ) ' '
  write ( *, * ) '  SEARCH RESULT:'
  if ( index > 0 ) then
    write ( *, * ) '    The value occurs in index ', index
  else
    write ( *, * ) '    The value does not occur in the array.'
  end if

  return
end
subroutine test109
!
!*******************************************************************************
!
!! TEST109 tests RVEC_ORDER_TYPE.
!
  integer, parameter :: n = 4
  integer, parameter :: ntest = 6
!
  integer itest
  integer j
  integer order
  real x(n,ntest)
!
  x(1,1) = 1.0E+00
  x(2,1) = 3.0E+00
  x(3,1) = 2.0E+00
  x(4,1) = 4.0E+00

  x(1,2) = 2.0E+00
  x(2,2) = 2.0E+00
  x(3,2) = 2.0E+00
  x(4,2) = 2.0E+00

  x(1,3) = 1.0E+00
  x(2,3) = 2.0E+00
  x(3,3) = 2.0E+00
  x(4,3) = 4.0E+00

  x(1,4) = 1.0E+00
  x(2,4) = 2.0E+00
  x(3,4) = 3.0E+00
  x(4,4) = 4.0E+00

  x(1,5) = 4.0E+00
  x(2,5) = 4.0E+00
  x(3,5) = 3.0E+00
  x(4,5) = 1.0E+00

  x(1,6) = 9.0E+00
  x(2,6) = 7.0E+00
  x(3,6) = 3.0E+00
  x(4,6) = 0.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST109'
  write ( *, * ) '  RVEC_ORDER_TYPE classifies a real vector as'
  write ( *, * ) '  -1: no order'
  write ( *, * ) '   0: all equal;'
  write ( *, * ) '   1: ascending;'
  write ( *, * ) '   2: strictly ascending;'
  write ( *, * ) '   3: descending;'
  write ( *, * ) '   4: strictly descending.'
  write ( *, * ) ' '

  do itest = 1, ntest

    call rvec_order_type ( n, x(1,itest), order )

    write ( *, * ) ' '
    write ( *, * ) 'The following vector has order type ', order
    write ( *, * ) ' '
    do j = 1, n
      write ( *, '(i8,g14.6)' ) j, x(j,itest)
    end do

  end do
   
  return
end
subroutine test110
!
!*******************************************************************************
!
!! TEST110 tests RVEC_PERMUTE.
!
  integer, parameter :: n = 5
!
  integer i
  integer perm(n)
  real x(n)
!
  perm(1) = 2
  perm(2) = 4
  perm(3) = 5
  perm(4) = 1
  perm(5) = 3

  x(1) = 1.0E+00
  x(2) = 2.0E+00
  x(3) = 3.0E+00
  x(4) = 4.0E+00
  x(5) = 5.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST110'
  write ( *, * ) '  RVEC_PERMUTE permutes a real vector in place.'
  write ( *, * ) ' '
  write ( *, * ) 'I, Perm(I), X(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(2i8,g14.6)' ) i, perm(i), x(i)
  end do

  call rvec_permute ( n, x, perm )

  call rvec_print ( n, x, '  Permuted array:' )

  return
end
subroutine test111
!
!*******************************************************************************
!
!! TEST111 tests RVEC_ROTATE.
!
  integer, parameter :: n = 5
!
  real a(n)
  integer m
!
  m = 2

  a(1) = 1.0E+00
  a(2) = 2.0E+00
  a(3) = 3.0E+00
  a(4) = 4.0E+00
  a(5) = 5.0E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST111'
  write ( *, * ) '  RVEC_ROTATE rotates a real vector in place.'
  write ( *, * ) ' '
  write ( *, * ) '  Rotate entries ', m, ' places to the right.'

  call rvec_print ( n, a, '  Original array:' )

  call rvec_rotate ( n, a, m )

  call rvec_print ( n, a, '  Rotated array:' )

  return
end
subroutine test112
!
!*******************************************************************************
!
!! TEST112 tests RVEC_REVERSE.
!
  integer, parameter :: n = 5
!
  real a(n)
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST112'
  write ( *, * ) '  RVEC_REVERSE reverses a real vector.'
  write ( *, * ) ' '
 
  call rvec_identity ( n, a )
 
  call rvec_print ( n, a, '  Original array:' )

  call rvec_reverse ( n, a )

  call rvec_print ( n, a, '  Reversed array:' )
 
  return
end
subroutine test113
!
!*******************************************************************************
!
!! TEST113 tests RVEC_SORT_HEAP_A.
!
  integer, parameter :: n = 20
!
  real a(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST113'
  write ( *, * ) '  RVEC_SORT_HEAP_A sorts a real vector.'
  write ( *, * ) ' '
 
  call rvec_random ( 0.0E+00, 3.0E+00 * real ( n ), n, a )
 
  call rvec_print ( n, a, '  Original array:' )

  call rvec_sort_heap_a ( n, a )

  call rvec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test114
!
!*******************************************************************************
!
!! TEST114 tests RVEC_SORT_HEAP_A;
!! TEST114 tests RVEC_UNIQ;
!! TEST114 tests RVEC_UNIQ_COUNT;
!! TEST114 tests RVEC_UNIQ_HIST.
!
  integer, parameter :: maxuniq = 30
  integer, parameter :: n = 30
!
  real a(n)
  real a_save(n)
  integer acount(maxuniq)
  real auniq(maxuniq)
  integer i
  integer nuniq
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST114'
  write ( *, * ) '  RVEC_SORT_HEAP_A sorts a real array;'
  write ( *, * ) '  RVEC_UNIQ finds the unique entries in'
  write ( *, * ) '    a sorted real array;'
  write ( *, * ) '  RVEC_UNIQ_COUNT counts the unique entries'
  write ( *, * ) '    of a sorted real array;'
  write ( *, * ) '  RVEC_UNIQ_HIST makes a histogram.'
  write ( *, * ) ' '

  call rvec_random ( 1.0E+00, real ( n ), n, a )

  a(1:n) = real ( int ( a(1:n) ) )

  a_save(1:n) = a(1:n)

  call rvec_print ( n, a, '  Unsorted array:' )
!
!  Test RVEC_UNIQ.
!
  call rvec_sort_heap_a ( n, a )
 
  call rvec_uniq ( n, a, nuniq )
 
  call rvec_print ( nuniq, a, '  Unique entries' )
!
!  Test RVEC_UNIQ_COUNT
!
  a(1:n) = a_save(1:n)
 
  call rvec_uniq_count ( n, a, nuniq )
 
  write ( *, * ) ' '
  write ( *, * ) '  RVEC_UNIQ_COUNT counts ', nuniq, ' unique entries in A.'
!
!  Test RVEC_UNIQ3
!
  a(1:n) = a_save(1:n)
 
  call rvec_sort_heap_a ( n, a )

  call rvec_uniq_hist ( n, a, maxuniq, nuniq, auniq, acount )
 
  write ( *, * ) ' '
  write ( *, * ) '  RVEC_UNIQ3 counts ', nuniq, ' unique entries.'
  write ( *, * ) ' '
  write ( *, * ) '  Value  Multiplicity'
  write ( *, * ) ' '
  do i = 1, nuniq
    write ( *, '(i8,g14.6,i8)' ) i, auniq(i), acount(i)
  end do

  return
end
subroutine test115
!
!*******************************************************************************
!
!! TEST115 tests RVEC_SORT_HEAP_INDEX_A.
!! TEST115 tests RVEC_SORT_HEAP_INDEX_D.
!
  integer, parameter :: n = 20
!
  real a(n)
  integer i
  integer indx(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST115'
  write ( *, * ) '  RVEC_SORT_HEAP_INDEX_A creates an ascending'
  write ( *, * ) '    sort index for a real array.'
  write ( *, * ) '  RVEC_SORT_HEAP_INDEX_D creates a descending'
  write ( *, * ) '    sort index for a real array.'

  call rvec_random ( 0.0E+00, 3.0E+00 * real ( n ), n, a )
 
  call rvec_print ( n, a, '  Unsorted array:' )

  call rvec_sort_heap_index_a ( n, a, indx )

  write ( *, * ) ' '
  write ( *, * ) '  After indexed ascending sort:'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), A(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(2i8,g14.6)' ) i, indx(i), a(i)
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Now use the index array to carry out the'
  write ( *, * ) '  permutation implicitly.'
  write ( *, * ) ' '
  write ( *, * ) '  INDX(I), A(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) indx(i), a(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Call RVEC_PERMUTE to carry out the permutation'
  write ( *, * ) '  explicitly.'
  write ( *, * ) ' '

  call rvec_permute ( n, a, indx )

  call rvec_print ( n, a, '  I, A(I)' )

  call rvec_sort_heap_index_d ( n, a, indx )

  write ( *, * ) ' '
  write ( *, * ) '  After indexed descending sort:'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), A(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(2i8,g14.6)' ) i, indx(i), a(i)
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Now use the index array to carry out the'
  write ( *, * ) '  permutation implicitly.'
  write ( *, * ) ' '
  write ( *, * ) '  INDX(I), ARRAY(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) indx(i), a(indx(i))
  end do

  return
end
subroutine test1154
!
!*******************************************************************************
!
!! TEST1154 tests RVEC_SORT_INSERT_A.
!
  integer, parameter :: n = 20
!
  real a(n)
  integer i
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1154'
  write ( *, * ) '  RVEC_SORT_INSERT_A ascending sorts a real array.'

  call rvec_random ( 0.0E+00, 3.0E+00 * real ( n ), n, a )

  call rvec_print ( n, a, '  Unsorted array:' )

  call rvec_sort_insert_a ( n, a )

  call rvec_print ( n, a, '  Sorted array:' )

  return
end
subroutine test1155
!
!*******************************************************************************
!
!! TEST1155 tests RVEC_SORT_INSERT_INDEX_A.
!
  integer, parameter :: n = 20
!
  real a(n)
  integer i
  integer indx(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1155'
  write ( *, * ) '  RVEC_SORT_INSERT_INDEX_A creates an ascending'
  write ( *, * ) '    sort index for a real array.'

  call rvec_random ( 0.0E+00, 3.0E+00 * real ( n ), n, a )

  call rvec_print ( n, a, '  Unsorted array:' )

  call rvec_sort_insert_index_a ( n, a, indx )

  write ( *, * ) ' '
  write ( *, * ) '  After indexed ascending sort:'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), A(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(2i8,g14.6)' ) i, indx(i), a(i)
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Now use the index array to carry out the'
  write ( *, * ) '  permutation implicitly.'
  write ( *, * ) ' '
  write ( *, * ) '  INDX(I), A(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i8,g14.6)' ) indx(i), a(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Call RVEC_PERMUTE to carry out the permutation'
  write ( *, * ) '  explicitly.'
  write ( *, * ) ' '

  call rvec_permute ( n, a, indx )

  call rvec_print ( n, a, '  I, A(I)' )

  return
end
subroutine test116
!
!*******************************************************************************
!
!! TEST116 tests RVEC_SPLIT_SORT.
!! TEST116 tests RVEC_SPLIT_UNSORT.
!
  integer, parameter :: n = 25
!
  real a(n)
  integer i
  integer i_gt
  integer i_lt
  integer isplit
  real split
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST116'
  write ( *, * ) '  RVEC_SPLIT_SORT splits a sorted vector into'
  write ( *, * ) '  entries less than and greater than a'
  write ( *, * ) '  splitting value.'
  write ( *, * ) '  RVEC_SPLIT_UNSORT splits an unsorted vector'
  write ( *, * ) '  in the same way.'
  write ( *, * ) ' '

  call rvec_random ( 0.0E+00, 10.0E+00, n, a )

  a(1:n) = real ( nint ( a(1:n) ) ) / 2.0E+00

  call rvec_sort_heap_a ( n, a )

  split = 0.5E+00 * ( a(1) + a(n) )

  call rvec_print ( n, a, '  The sorted array:' )

  write ( *, * ) ' '
  write ( *, * ) '  Splitting value is ', split
  write ( *, * ) ' '

  call rvec_split_sort ( n, a, split, i_lt, i_gt )

  write ( *, * ) '  Lower index I_LT = ', i_lt
  write ( *, * ) '  Upper index I_GT = ', i_gt

  write ( *, * ) ' '
  write ( *, * ) '  Now repeat test with RVEC_SPLIT_UNSORT.'
  write ( *, * ) ' '
  call rvec_permute_random ( n, a )
 
  call rvec_print ( n, a, '  The shuffled array:' )
 
  call rvec_split_unsort ( n, a, split, isplit )
 
  call rvec_print ( n, a, '  The split array:' )

  write ( *, * ) ' '
  write ( *, * ) '  Array entries <= SPLIT up to index ', isplit
 
  return
end
subroutine test1165
!
!*******************************************************************************
!
!! TEST1165 tests RVEC2_SORT_A.
!! TEST1165 tests RVEC2_SORT_D.
!! TEST1165 tests RVEC2_UNIQ.
!
  integer, parameter :: n = 10
!
  real a1(n)
  real a2(n)
  integer i
  integer nuniq
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1165'
  write ( *, * ) '  For a pair of real vectors:'
  write ( *, * ) '  RVEC2_SORT_A ascending sorts;'
  write ( *, * ) '  RVEC2_SORT_D descending sorts;'
  write ( *, * ) '  RVEC2_UNIQ counts unique entries.'
  write ( *, * ) ' '

  call rvec_random ( 1.0E+00, 3.0E+00, n, a1 )

  call rvec_random ( 5.0E+00, 10.0E+00, n, a2 )

  a1(3) = a1(1)
  a2(3) = a2(1)

  a1(6) = a1(2)
  a2(6) = a2(2)

  a1(9) = a1(1)
  a2(9) = a2(1)

  call rvec2_print ( n, a1, a2, '  The pair of arrays:' )

  call rvec2_sort_a ( n, a1, a2 )

  call rvec2_print ( n, a1, a2, '  Arrays after ascending sort:' )

  call rvec2_sort_d ( n, a1, a2 )

  call rvec2_print ( n, a1, a2, '  Arrays after descending sort:' )

  call rvec2_uniq ( n, a1, a2, nuniq )

  call rvec2_print ( nuniq, a1, a2, '  UNIQed array:' )

  return
end
subroutine test126
!
!*******************************************************************************
!
!! TEST126 tests RVEC2_SORT_HEAP_INDEX_A.
!
  integer, parameter :: n = 20
!
  integer i
  integer indx(n)
  real temp
  real x(n)
  real y(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST126'
  write ( *, * ) '  RVEC2_SORT_HEAP_INDEX_A creates a sort index'
  write ( *, * ) '  for an (X,Y) array.'
  write ( *, * ) ' '
 
  do i = 1, n

    call random_number ( harvest = temp )
    x(i) = real ( int ( real ( n ) * temp ) ) / real ( n )

    call random_number ( harvest = temp )
    y(i) = real ( int ( real ( n ) * temp ) ) / real ( n )

  end do
 
  write ( *, * ) '  The unsorted array:'
  write ( *, * ) ' '
  write ( *, * ) '  I, X(I), Y(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,6x,2g14.6)' ) i, x(i), y(i)
  end do

  call rvec2_sort_heap_index_a ( n, x, y, indx )

  write ( *, * ) ' '
  write ( *, * ) '  After sorting:'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), X(I), Y(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(2i6,2g14.6)' ) i, indx(i), x(i), y(i)
  end do
 
  write ( *, * ) ' '
  write ( *, * ) '  Now use the index array to carry out the'
  write ( *, * ) '  permutation implicitly.'
  write ( *, * ) ' '
  write ( *, * ) '  I, INDX(I), X(INDX(I)), Y(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(2i6,2g14.6)' ) i, indx(i), x(indx(i)), y(indx(i))
  end do

  write ( *, * ) ' '
  write ( *, * ) '  RVEC_PERMUTE carries out the permutation.'

  call rvec_permute ( n, x, indx )
  call rvec_permute ( n, y, indx )

  write ( *, * ) ' '
  write ( *, * ) '  I, X(I), Y(I)'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i6,6x,2g14.6)' ) i, x(i), y(i)
  end do

  return
end
subroutine test117
!
!*******************************************************************************
!
!! TEST117 tests RVEC2_SUM_IMAX.
!
  integer, parameter :: n = 10
!
  real a(n)
  real b(n)
  integer ival
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST117'
  write ( *, * ) '  For 2 real vectors:'
  write ( *, * ) '  RVEC2_SUM_IMAX: index of the sum vector'
  write ( *, * ) '    with maximum value.'
  write ( *, * ) ' '

  call rvec_random ( 0.0E+00, 10.0E+00, n, a )
  call rvec_random ( 0.0E+00,  5.0E+00, n, b )

  call rvec2_print ( n, a, b, '  The pair of vectors:' )

  write ( *, * ) ' '

  call rvec2_sum_imax ( n, a, b, ival )

  write ( *, * ) '  Index of maximum in A+B: ', ival

  return
end
subroutine test1174
!
!*******************************************************************************
!
!! TEST1174 tests RVEC3_INDEX_SEARCH.
!! TEST1174 tests RVEC3_INDEX_INSERT_UNIQUE.
!
  integer, parameter :: maxn = 30
!
  integer equal
  integer i
  integer ierror
  integer indx(maxn)
  integer ival
  integer j
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
  n = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST1174'
  write ( *, * ) '  RVEC3_INDEX_INSERT_UNIQUE inserts unique values into an'
  write ( *, * ) '    index sorted array.'
  write ( *, * ) '  RVEC3_INDEX_SEARCH searches for an entry with a given value.'
  write ( *, * ) ' '
  write ( *, * ) '  Generate ', maxn, ' random values:'
  write ( *, * ) ' '
  do i = 1, maxn
    call r_random ( 1.0E+00, 4.0E+00, xval )
    xval = real ( nint ( xval ) )
    call r_random ( 1.0E+00, 3.0E+00, yval )
    yval = real ( nint ( yval ) )
    call r_random ( 1.0E+00, 4.0E+00, zval )
    zval = real ( nint ( zval ) )
    write ( *, '(4x,3f6.2)' ) xval, yval, zval
    call rvec3_index_insert_unique ( maxn, n, x, y, z, indx, &
      xval, yval, zval, ival, ierror )
    write ( *, * ) 'IVAL = ', ival
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Indexed list of unique entries:'
  write ( *, * ) ' '
  write ( *, '(a)' ) '  I  INDX(I)  X(I)  X(INDX(I))'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(i3,6x,i3,3x,f6.2,9x,f6.2,9x,f6.2)' ) i, indx(i), x(indx(i)), &
      y(indx(i)), z(indx(i))
  end do

  return
end
subroutine test1175
!
!*******************************************************************************
!
!! TEST1175 tests SGE_HESS.
!
  integer, parameter :: n = 3
  integer, parameter :: lda = n
!
  real h(lda,n)
  real x(n)
!
  external test1175_f
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1175'
  write ( *, * ) '  SGE_HESS estimates the Hessian matrix'
  write ( *, * ) '    of a scalar function.'

  x(1) = 1.0E+00
  x(2) = 2.0E+00
  x(3) = 3.0E+00

  call sge_hess ( test1175_f, lda, n, x, h )

  call rmat_print ( lda, n, n, h, '  Estimated jacobian:' )

  call test1175_hess ( lda, n, x, h )

  call rmat_print ( lda, n, n, h, '  Exact jacobian:' )

  return
end
subroutine test1175_f ( n, x, f )
!
!*******************************************************************************
!
!! TEST1175_F is a sample nonlinear function for treatment by SGE_JAC.
!
  integer n
!
  real f
  real x(n)
!
  f = x(1)**2 + x(1) * x(2) + x(2) * cos ( 10.0 * x(3) )

  return
end
subroutine test1175_hess ( lda, n, x, h )
!
!*******************************************************************************
!
!! TEST1175_HESS is the exact Hessian of TEST1175_F.
!
  integer lda
  integer n
!
  real h(lda,n)
  real x(n)
!
  h(1,1) = 2.0E+00
  h(1,2) = 1.0E+00
  h(1,3) = 0.0E+00

  h(2,1) = 1.0E+00
  h(2,2) = 0.0E+00
  h(2,3) = - 10.0E+00 * sin ( 10.0E+00 * x(3) )

  h(3,1) = 0.0E+00
  h(3,2) = - 10.0E+00 * sin ( 10.0E+00 * x(3) )
  h(3,3) = - 100.0E+00 * x(2) * cos ( 10.0E+00 * x(3) )

  return
end
subroutine test118
!
!*******************************************************************************
!
!! TEST118 tests SGE_JAC.
!
  integer, parameter :: m = 3
  integer, parameter :: n = 4
  integer, parameter :: lda = m
!
  real eps
  real fprime(lda,n)
  real x(n)
!
  external test118_f
!
  eps = 0.00001E+00

  write ( *, * ) ' '
  write ( *, * ) 'TEST118'
  write ( *, * ) '  SGE_JAC estimates the M by N jacobian matrix'
  write ( *, * ) '    of a nonlinear function.'

  x(1) = 1.0E+00
  x(2) = 2.0E+00
  x(3) = 3.0E+00
  x(4) = 4.0E+00

  call sge_jac ( eps, fprime, test118_f, lda, m, n, x  )

  call rmat_print ( lda, m, n, fprime, '  Estimated jacobian:' )

  call test118_jac ( lda, m, n, fprime, x )

  call rmat_print ( lda, m, n, fprime, '  Exact jacobian:' )

  return
end
subroutine test118_f ( m, n, x, f )
!
!*******************************************************************************
!
!! TEST118_F is a sample nonlinear function for treatment by SGE_JAC.
!
  integer m
  integer n
!
  real f(m)
  real x(n)
!
  f(1) = sin ( x(1) * x(2) )
  f(2) = sqrt ( 1.0E+00 + x(1)**2 ) + x(3)
  f(3) = x(1) + 2.0E+00 * x(2) + 3.0E+00 * x(3) + 4.0E+00 * x(4)

  return
end
subroutine test118_jac ( lda, m, n, fprime, x )
!
!*******************************************************************************
!
!! TEST118_JAC is the exact jacobian of TEST118_F.
!
  integer lda
  integer m
  integer n
!
  real fprime(lda,n)
  real x(n)
!
  if ( lda < m ) then
    write ( *, * ) ' '
    write ( *, * ) 'FP_JAC - Fatal error!'
    write ( *, * ) '  LDA < M.'
    stop
  end if

  fprime(1,1) = cos ( x(1) * x(2) ) * x(2)
  fprime(1,2) = cos ( x(1) * x(2) ) * x(1)
  fprime(1,3) = 0.0E+00
  fprime(1,4) = 0.0E+00

  fprime(2,1) = x(1) / sqrt ( 1.0E+00 + x(1)**2 )
  fprime(2,2) = 0.0E+00
  fprime(2,3) = 1.0E+00
  fprime(2,4) = 0.0E+00

  fprime(3,1) = 1.0E+00
  fprime(3,2) = 2.0E+00
  fprime(3,3) = 3.0E+00
  fprime(3,4) = 4.0E+00

  return
end
subroutine test119
!
!*******************************************************************************
!
!! TEST119 tests SGE_SOLVE.
!
  integer, parameter :: lda = 4
!
  real a(lda,lda)
  real b(4)
  integer ierror
  integer itest
  integer n
  real x(4)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST119'
  write ( *, * ) '  SGE_SOLVE is a linear solver.'
  write ( *, * ) ' '
!
  do itest = 1, 4
 
    if ( itest == 1 ) then

      n = 2
 
      a(1,1) = 1.0E+00
      a(1,2) = 2.0E+00
      a(2,1) = 3.0E+00
      a(2,2) = 4.0E+00
     
      b(1) = 5.0E+00
      b(2) = 11.0E+00
 
    else if ( itest == 2 ) then
 
      n = 3
     
      a(1,1) = 2.0E+00
      a(1,2) = 1.0E+00
      a(1,3) = 1.0E+00
     
      a(2,1) = 1.0E+00
      a(2,2) = 1.0E+00
      a(2,3) = 0.0E+00
     
      a(3,1) = 1.0E+00
      a(3,2) = 0.0E+00
      a(3,3) = 1.0E+00
     
      b(1) = 4.0E+00
      b(2) = 2.0E+00
      b(3) = 2.0E+00
 
    else if ( itest == 3 ) then
 
      n = 4
     
      a(1,1) = 1.0E+00
      a(1,2) = 0.0E+00
      a(1,3) = 0.0E+00
      a(1,4) = 1.0E+00
     
      a(2,1) = 2.0E+00
      a(2,2) = 1.0E+00
      a(2,3) = 0.0E+00
      a(2,4) = 3.0E+00
     
      a(3,1) = 1.0E+00
      a(3,2) = 2.0E+00
      a(3,3) = 3.0E+00
      a(3,4) = 0.0E+00
     
      a(4,1) = 3.0E+00
      a(4,2) = 1.0E+00
      a(4,3) = 2.0E+00
      a(4,4) = 1.0E+00
     
      b(1) = 5.0E+00
      b(2) = 16.0E+00
      b(3) = 14.0E+00
      b(4) = 15.0E+00
 
    else if ( itest == 4 ) then
 
      n = 3
     
      a(1,1) = 2.0E+00
      a(1,2) = 4.0E+00
      a(1,3) = 1.0E+00
     
      a(2,1) = 1.0E+00
      a(2,2) = 2.0E+00
      a(2,3) = 4.0E+00
     
      a(3,1) = 3.0E+00
      a(3,2) = 6.0E+00
      a(3,3) = 5.0E+00
     
      b(1) = 13.0E+00
      b(2) = 17.0E+00
      b(3) = 20.0E+00
     
    end if
 
    call rvec_print ( n, b, '  Right hand side:' )
 
    call sge_solve ( a, b, lda, n, x, ierror )
 
    write ( *, * ) ' '
    if ( ierror == 0 ) then
      write ( *, * ) '  The system is nonsingular.'
    else if ( ierror == 1 ) then
      write ( *, * ) '  The system is singular, but consistent.'
    else if ( ierror == 2 ) then
      write ( *, * ) '  The system is singular and inconsistent.'
    end if
   
    call rvec_print ( n, x, '  Computed solution:' )
 
  end do
 
  return
end
subroutine test120
!
!*******************************************************************************
!
!! TEST120 tests SORT_HEAP_EXTERNAL.
!
  integer,  parameter :: n = 20
!
  integer a(n)
  integer i
  integer indx
  integer isgn
  integer j
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST120'
  write ( *, * ) '  SORT_HEAP_EXTERNAL sorts objects externally.'
  write ( *, * ) ' '

  indx = 0
  i = 0
  j = 0
  isgn = 0
 
  call ivec_random ( n, 1, n, a )
 
  call ivec_print ( n, a, '  Unorted array:' )
 
  do

    call sort_heap_external ( n, indx, i, j, isgn )
 
    if ( indx < 0 ) then

      isgn = 1
      if ( a(i) <= a(j) ) then
        isgn = -1
      end if

    else if ( indx > 0 ) then

      call i_swap ( a(i), a(j) )

    else

      exit

    end if

  end do

  call ivec_print ( n, a, '  Sorted array:' )
 
  return
end
subroutine test121
!
!*******************************************************************************
!
!! TEST121 tests SVEC_REVERSE.
!! TEST121 tests SVEC_SORT_HEAP_A.
!
  integer, parameter :: n = 10
!
  character ( len = 10 ) carray(n)
  integer i
!
  carray(1) = 'FRED'
  carray(2) = 'fred'
  carray(3) = 'Abacus'
  carray(4) = 'beetles'
  carray(5) = 'XYLOPHONE'
  carray(6) = 'banana'
  carray(7) = 'goofball'
  carray(8) = 'abbot'
  carray(9) = 'BARBECUE'
  carray(10) = 'abbots'
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST121'
  write ( *, * ) '  SVEC_SORT_HEAP_A sorts a string vector.'
  write ( *, * ) '  SVEC_REVERSE reverses a string vector.'
  write ( *, * ) ' '
  write ( *, * ) '  Unsorted list:'
  write ( *, * ) ' '
 
  do i = 1, n
    write ( *, '(5x,a)' ) carray(i)
  end do
 
  call svec_sort_heap_a ( n, carray )
 
  write ( *, * ) ' '
  write ( *, * ) '  Sorted list:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(5x,a)' ) carray(i)
  end do
 
  call svec_reverse ( n, carray )
 
  write ( *, * ) ' '
  write ( *, * ) '  Reversed sorted list:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(5x,a)' ) carray(i)
  end do

  return
end
subroutine test122
!
!*******************************************************************************
!
!! TEST122 tests SVEC_SORT_HEAP_A;
!! TEST122 tests SVEC_MERGE_A.
!! TEST122 tests SVEC_SEARCH_BINARY_A.
!
  integer, parameter :: na = 10
  integer, parameter :: nb = 10
!
  character ( len = 4 ) a(na)
  character ( len = 4 ) b(nb)
  character ( len = 4 ) c(na+nb)
  integer i
  integer indx
  integer j
  integer nc
  character ( len = 4 ) string
  real u
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST122'
  write ( *, * ) '  For ascending order:'
  write ( *, * ) '  SVEC_SORT_HEAP_A sorts a character array;'
  write ( *, * ) '  SVEC_MERGE_A merges two sorted character '
  write ( *, * ) '    arrays into a single sorted array.'
  write ( *, * ) '  SVEC_SEARCH_BINARY_A searches a string array for'
  write ( *, * ) '    a particular value.'
  write ( *, * ) ' '
 
  do i = 1, na
    do j = 1, 4
      call ch_random ( 'A', 'E', a(i) (j:j) )
    end do
  end do
 
  call svec_sort_heap_a ( na, a )
 
  write ( *, * ) ' '
  write ( *, * ) '  Sorted vector A:'
  write ( *, * ) ' '
 
  do i = 1, na
    write ( *, * ) a(i)
  end do

  do i = 1, nb
    do j = 1, 4
      call ch_random ( 'B', 'F', b(i) (j:j) )
    end do
  end do
 
  call svec_sort_heap_a ( nb, b )
 
  write ( *, * ) ' '
  write ( *, * ) '  Sorted vector B:'
  write ( *, * ) ' '
 
  do i = 1, nb
    write ( *, * ) b(i)
  end do
 
  call svec_merge_a ( na, a, nb, b, nc, c )

  write ( *, * ) ' '
  write ( *, * ) '  Merged output vector C = A + B:'
  write ( *, * ) ' '
 
  do i = 1, nc
    write ( *, * ) c(i)
  end do

  string = a(2)

  write ( *, * ) ' '
  write ( *, * ) '  Search C for value ' // string
  write ( *, * ) ' '

  call svec_search_binary_a ( nc, c, string, indx )

  if ( indx == 0 ) then
    write ( *, * ) '  The value does not occur'
  else
    write ( *, * ) '  The value occurs at index ', indx
  end if
 
  return
end
subroutine test123
!
!*******************************************************************************
!
!! TEST123 tests SVEC_SORT_HEAP_A;
!! TEST123 tests SVEC_UNIQ.
!
  integer, parameter :: n = 10
!
  character ( len = 3 ) a(n)
  integer i
  integer nuniq
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST123'
  write ( *, * ) '  SVEC_SORT_HEAP_A sorts a character array;'
  write ( *, * ) '  SVEC_UNIQ finds the unique entries in a'
  write ( *, * ) '    sorted character array.'
  write ( *, * ) ' '
 
  a(1) = 'Cat'
  a(2) = 'Bat'
  a(3) = 'Mat'
  a(4) = 'Tab'
  a(5) = 'Ax'
  a(6) = 'Ax'
  a(7) = 'Tab'
  a(8) = 'Pyx'
  a(9) = 'Ax'
  a(10) = 'Bat'
 
  call svec_sort_heap_a ( n, a )
 
  write ( *, * ) ' '
  write ( *, * ) '  Input vector A:'
  write ( *, * ) ' '
 
  do i = 1, n
    write ( *, '(5x,a)' ) a(i)
  end do
 
  call svec_uniq ( n, a, nuniq )
 
  write ( *, * ) ' '
  write ( *, * ) '  Unique entries:'
  write ( *, * ) ' '
 
  do i = 1, nuniq
    write ( *, '(5x,a)' ) a(i)
  end do
 
  return
end
subroutine test124
!
!*******************************************************************************
!
!! TEST124 tests SVECI_SEARCH_BINARY_A.
!! TEST124 tests SVECI_SORT_HEAP_A.
!
  integer, parameter :: n = 10
!
  character ( len = 10 ) carray(n)
  integer i
  integer indx
  character ( len = 10 ) string
!
  carray(1) = 'FRED'
  carray(2) = 'fred'
  carray(3) = 'Abacus'
  carray(4) = 'beetles'
  carray(5) = 'XYLOPHONE'
  carray(6) = 'banana'
  carray(7) = 'goofball'
  carray(8) = 'abbot'
  carray(9) = 'BARBECUE'
  carray(10) = 'abbots'
 
  write ( *, * ) ' '
  write ( *, * ) 'TEST124'
  write ( *, * ) '  For implicitly capitalized strings,'
  write ( *, * ) '  SVECI_SORT_HEAP_A sorts;'
  write ( *, * ) '  SVECI_SEARCH_BINARY_A searches.'
  write ( *, * ) ' '
  write ( *, * ) '  Unsorted list:'
  write ( *, * ) ' '
 
  do i = 1, n
    write ( *, '(5x,a)' ) carray(i)
  end do
 
  call sveci_sort_heap_a ( n, carray )
 
  write ( *, * ) ' '
  write ( *, * ) '  Sorted list:'
  write ( *, * ) ' '
  do i = 1, n
    write ( *, '(5x,a)' ) carray(i)
  end do
 
  string = 'ABBoT'

  write ( *, * ) ' '
  write ( *, * ) '  Now search for the string ' // string

  call sveci_search_binary_a ( n, carray, string, indx )
 
  write ( *, * ) ' '
  if ( indx == 0 ) then
    write ( *, * ) '  The search string does not occur.'
  else
    write ( *, * ) '  The search string occurs in index ', indx
  end if

  return
end
subroutine test1245
!
!*******************************************************************************
!
!! TEST1245 tests UNIFORM_01_SAMPLE
!
  integer i
  integer iseed
  real mean
  integer, parameter :: n = 1000
  real uniform_01_sample
  real variance
  real x(n)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST1245'
  write ( *, * ) '  UNIFORM_01_SAMPLE samples a uniform random'
  write ( *, * ) '  distribution in [0,1].'

  call get_seed ( iseed )

  write ( *, * ) ' '
  write ( *, * ) '  Starting with seed = ', iseed

  do i = 1, n
    x(i) = uniform_01_sample ( iseed )
  end do

  call rvec_mean ( n, x, mean )

  call rvec_variance ( n, x, variance )

  write ( *, * ) '  Number of values computed was N = ', n
  write ( *, * ) '  Average value was ', mean
  write ( *, * ) '  Variance was ', variance

  return
end
subroutine test127
!
!*******************************************************************************
!
!! TEST127 tests C_CUBE_ROOT.
!
  integer, parameter :: ntest = 5
!
  integer i
  complex x(ntest)
  complex y
  complex z
  complex c_cube_root
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST127'
  write ( *, * ) '  C_CUBE_ROOT finds complex cube roots.'
 
  x(1) = cmplx (  0.0E+00,  0.0E+00 )
  x(2) = cmplx ( -1.0E+00,  0.0E+00 )
  x(3) = cmplx (  0.0E+00,  8.0E+00 )
  x(4) = cmplx (  1.0E+00,  1.0E+00 ) / sqrt ( 2.0E+00 )
  x(5) = cmplx (  0.0E+00, -1.0E+00 )

  write ( *, * ) ' '
  write ( *, * ) '          X                  Y=C_CUBE_ROOT(X)            Y**3'
  write ( *, * ) ' '

  do i = 1, ntest
 
    y = c_cube_root ( x(i) )
    z = y**3

    write ( *, '(3(f10.5,f10.5,5x))' ) x(i), y, z
 
  end do

  return
end
subroutine test128
!
!*******************************************************************************
!
!! TEST128 tests R_ZETA.
!
  integer i
  real p
  real r_zeta
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST128'
  write ( *, * ) '  R_ZETA computes the Riemann Zeta function.'
  write ( *, * ) ' '
  write ( *, * ) '  P    Zeta(P)'
  write ( *, * ) ' '
  do i = 3, 12
    p = real ( i ) / 2.0E+00
    write ( *, '(2g14.6)' ) p, r_zeta(p)
  end do

  return
end

