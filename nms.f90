function alngam ( x )
!
!*******************************************************************************
!
!! ALNGAM computes the log of the absolute value of the gamma function.
!
!
!  Definition:
!
!    The Gamma function is defined as
!
!      GAMMA(Z) = INTEGRAL ( 0 <= T < Infinity) T**(Z-1) EXP(-T) DT
!
!    If Z is a positive integer, GAMMA(Z) = (Z-1)!, the factorial.
!
!    There is a special value:
!
!      GAMMA(0.5) = SQRT(PI).
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Modified:
!
!    31 May 2000
!
!  Parameters:
!
!    Input, real X, the argument of the gamma function.
!
!    Output, real ALNGAM, the logarithm of the absolute value of GAMMA(X).
!
  real alngam
  real, save :: dxrel = 0.0E+00
  real gamma
  real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510E+00
  real r1mach
  real r9lgmc
  real sinpiy
  real, parameter :: sq2pil = 0.91893853320467274E+00
  real, parameter :: sqpi2l = 0.22579135264472743E+00
  real x
  real, save :: xmax = 0.0E+00
  real y
!
  if ( xmax == 0.0E+00 ) then
    xmax = huge ( xmax ) / log ( huge ( xmax ) )
    dxrel = sqrt ( epsilon ( dxrel ) )
  end if

  y = abs ( x )

  if ( y <= 10.0E+00 ) then

    alngam = log ( abs ( gamma ( x ) ) )
    return

  end if

  if ( y > xmax ) then
    write ( *, * ) ' '
    write ( *, * ) 'ALNGAM - Fatal error!'
    write ( *, * ) '  |X| is so big that ALNGAM will overflow.'
    stop
  end if

  if ( x > 0.0E+00 ) then
    alngam = sq2pil + ( x - 0.5E+00 )* log ( x ) - x + r9lgmc ( y )
    return
  end if

  sinpiy = abs ( sin ( pi * y ) )

  if ( sinpiy == 0.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ALNGAM - Fatal error!'
    write ( *, * ) '  X is a negative integer.'
    stop
  end if

  if ( abs ( ( x - aint ( x - 0.5E+00 ) ) / x ) < dxrel ) then
    write ( *, * ) ' '
    write ( *, * ) 'ALNGAM - Warning:'
    write ( *, * ) '  Answer has less than half usual precision.'
    write ( *, * ) '  X is very near a negative integer.'
  end if

  alngam = sqpi2l + ( x - 0.5E+00 ) * log ( y ) - x - log ( sinpiy ) &
    - r9lgmc ( y )

  return
end
subroutine asyjy ( funjy, x, fnu, flgjy, in, y, wk, iflw )
!
!*******************************************************************************
!
!! ASYJY computes high order Bessel functions J and Y.
!
!
!  Description:
!
!    ASYJY implements the uniform asymptotic expansion of
!    the J and Y Bessel functions for FNU >= 35 and real
!    X > 0.0.  The forms are identical except for a change
!    in sign of some of the terms.  This change in sign is
!    accomplished by means of the flag FLGJY = 1 or -1. On
!    flgjy = 1 the airy functions ai(x) and dai(x) are
!    supplied by the external function jairy, and on
!    flgjy = -1 the airy functions bi(x) and dbi(x) are
!    supplied by the external funtion yairy.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Author:
!
!    D E Amos
!
!  Parameters:
!
!    Input, external FUNJY, is the function JAIRY or YAIRY.
!
!    Input, real X, the argument, which must be greater than 0.
!
!    Input, real FNU, the order of the first Bessel function.
!    FNU is generally at least 35.
!
!    Input, real FLGJY, a selection flag
!     1.0E+00 gives the j function
!    -1.0E+00 gives the y function
!
!    Input, integer IN, the number of functions desired, which should be
!    1 or 2.
!
!    Output, real Y(IN), contains the desired function values.
!
!    Output, integer IFLW, a flag indicating underflow or overflow
!    return variables for besj only.
!
!    Output, real WK(7), contains the following values:
!
!      wk(1) = 1 - (x/fnu)**2 = w**2
!      wk(2) = sqrt(abs(wk(1)))
!      wk(3) = abs(wk(2) - atan(wk(2)))  or
!              abs(ln((1 + wk(2))/(x/fnu)) - wk(2))
!            = abs((2/3)*zeta**(3/2))
!      wk(4) = fnu*wk(3)
!      wk(5) = (1.5*wk(3)*fnu)**(1/3) = sqrt(zeta)*fnu**(1/3)
!      wk(6) = sign(1.,w**2)*wk(5)**2 = sign(1.,w**2)*zeta*fnu**(2/3)
!      wk(7) = fnu**(1/3)
!
  real abw2
  real akm
  real alfa(26,4)
  real alfa1
  real alfa2
  real ap
  real ar(8)
  real asum
  real az
  real beta(26,5)
  real beta1
  real beta2
  real beta3
  real br(10)
  real bsum
  real c(65)
  real con1
  real con2
  real con548
  real cr(10)
  real crz32
  real dfi
  real elim
  real dr(10)
  real fi
  real flgjy
  real fn
  real fnu
  real fn2
  real gama(26)
  integer i
  integer i1mach
  integer iflw
  integer in
  integer j
  integer jn
  integer jr
  integer ju
  integer k
  integer kb
  integer klast
  integer kmax(5)
  integer kp1
  integer ks
  integer ksp1
  integer kstemp
  integer l
  integer lr
  integer lrp1
  real phi
  real rcz
  real rden
  real relb
  real rfn2
  real rtz
  real rzden
  real sa
  real sb
  real suma
  real sumb
  real s1
  real ta
  real tau
  real tb
  real tfn
  real tol
  real, save :: tols = -6.90775527898214E+00
  real t2
  real upol(10)
  real wk(*)
  real x
  real xx
  real y(*)
  real z
  real z32
  real r1mach

  dimension alfa1(26,2), alfa2(26,2)
  dimension beta1(26,2), beta2(26,2), beta3(26,1)
!
  external funjy
!
  equivalence (alfa(1,1),alfa1(1,1))
  equivalence (alfa(1,3),alfa2(1,1))
  equivalence (beta(1,1),beta1(1,1))
  equivalence (beta(1,3),beta2(1,1))
  equivalence (beta(1,5),beta3(1,1))
!
  data con1/6.66666666666667e-01/
  data con2/3.33333333333333e-01/
  data con548/1.04166666666667e-01/
  data  ar(1),  ar(2),  ar(3),  ar(4),  ar(5),  ar(6),  ar(7), ar(8) &
    / 8.35503472222222e-02, 1.28226574556327e-01, &
      2.91849026464140e-01, 8.81627267443758e-01, 3.32140828186277, &
      1.49957629868626e+01, 7.89230130115865e+01, 4.74451538868264e+02/

  data  br(1), br(2), br(3), br(4), br(5), br(6), br(7), br(8), br(9), br(10) &
    /-1.45833333333333e-01,-9.87413194444444e-02, &
     -1.43312053915895e-01,-3.17227202678414e-01,-9.42429147957120e-01, &
     -3.51120304082635,-1.57272636203680e+01,-8.22814390971859e+01, &
     -4.92355370523671e+02,-3.31621856854797e+03/

  data c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10), &
     c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18), &
     c(19), c(20), c(21), c(22), c(23), c(24)/ &
     -2.08333333333333e-01,        1.25000000000000e-01, &
      3.34201388888889e-01,       -4.01041666666667e-01, &
      7.03125000000000e-02,       -1.02581259645062, &
      1.84646267361111,       -8.91210937500000e-01, &
      7.32421875000000e-02,        4.66958442342625, &
     -1.12070026162230e+01,        8.78912353515625, &
     -2.36408691406250,        1.12152099609375e-01, &
     -2.82120725582002e+01,        8.46362176746007e+01, &
     -9.18182415432400e+01,        4.25349987453885e+01, &
     -7.36879435947963,        2.27108001708984e-01, &
      2.12570130039217e+02,       -7.65252468141182e+02, &
      1.05999045252800e+03,       -6.99579627376133e+02/

  data c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32), &
          c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40), &
          c(41), c(42), c(43), c(44), c(45), c(46), c(47), c(48)/ &
             2.18190511744212e+02,       -2.64914304869516e+01, &
             5.72501420974731e-01,       -1.91945766231841e+03, &
             8.06172218173731e+03,       -1.35865500064341e+04, &
             1.16553933368645e+04,       -5.30564697861340e+03, &
             1.20090291321635e+03,       -1.08090919788395e+02, &
             1.72772750258446,        2.02042913309661e+04, &
            -9.69805983886375e+04,        1.92547001232532e+05, &
            -2.03400177280416e+05,        1.22200464983017e+05, &
            -4.11926549688976e+04,        7.10951430248936e+03, &
            -4.93915304773088e+02,        6.07404200127348, &
            -2.42919187900551e+05,        1.31176361466298e+06, &
            -2.99801591853811e+06,        3.76327129765640e+06/

  data c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56), &
          c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64), &
          c(65)/ &
            -2.81356322658653e+06,        1.26836527332162e+06, &
            -3.31645172484564e+05,        4.52187689813627e+04, &
            -2.49983048181121e+03,        2.43805296995561e+01, &
             3.28446985307204e+06,       -1.97068191184322e+07, &
             5.09526024926646e+07,       -7.41051482115327e+07, &
             6.63445122747290e+07,       -3.75671766607634e+07, &
             1.32887671664218e+07,       -2.78561812808645e+06, &
             3.08186404612662e+05,       -1.38860897537170e+04, &
             1.10017140269247e+02/

  data alfa1(1,1), alfa1(2,1), alfa1(3,1), alfa1(4,1), alfa1(5,1), &
          alfa1(6,1), alfa1(7,1), alfa1(8,1), alfa1(9,1), alfa1(10,1), &
          alfa1(11,1),alfa1(12,1),alfa1(13,1),alfa1(14,1),alfa1(15,1), &
          alfa1(16,1),alfa1(17,1),alfa1(18,1),alfa1(19,1),alfa1(20,1), &
          alfa1(21,1),alfa1(22,1),alfa1(23,1),alfa1(24,1),alfa1(25,1), &
          alfa1(26,1)     /-4.44444444444444e-03,-9.22077922077922e-04, &
     -8.84892884892885e-05, 1.65927687832450e-04, 2.46691372741793e-04, &
      2.65995589346255e-04, 2.61824297061501e-04, 2.48730437344656e-04, &
      2.32721040083232e-04, 2.16362485712365e-04, 2.00738858762752e-04, &
      1.86267636637545e-04, 1.73060775917876e-04, 1.61091705929016e-04, &
      1.50274774160908e-04, 1.40503497391270e-04, 1.31668816545923e-04, &
      1.23667445598253e-04, 1.16405271474738e-04, 1.09798298372713e-04, &
      1.03772410422993e-04, 9.82626078369363e-05, 9.32120517249503e-05, &
      8.85710852478712e-05, 8.42963105715700e-05, 8.03497548407791e-05/

  data alfa1(1,2), alfa1(2,2), alfa1(3,2), alfa1(4,2), alfa1(5,2), &
          alfa1(6,2), alfa1(7,2), alfa1(8,2), alfa1(9,2), alfa1(10,2), &
          alfa1(11,2),alfa1(12,2),alfa1(13,2),alfa1(14,2),alfa1(15,2), &
          alfa1(16,2),alfa1(17,2),alfa1(18,2),alfa1(19,2),alfa1(20,2), &
          alfa1(21,2),alfa1(22,2),alfa1(23,2),alfa1(24,2),alfa1(25,2), &
          alfa1(26,2)     / 6.93735541354589e-04, 2.32241745182922e-04, &
     -1.41986273556691e-05,-1.16444931672049e-04,-1.50803558053049e-04,&
     -1.55121924918096e-04,-1.46809756646466e-04,-1.33815503867491e-04, &
     -1.19744975684254e-04,-1.06184319207974e-04,-9.37699549891194e-05, &
     -8.26923045588193e-05,-7.29374348155221e-05,-6.44042357721016e-05, &
     -5.69611566009369e-05,-5.04731044303562e-05,-4.48134868008883e-05, &
     -3.98688727717599e-05,-3.55400532972042e-05,-3.17414256609022e-05, &
     -2.83996793904175e-05,-2.54522720634871e-05,-2.28459297164725e-05, &
     -2.05352753106481e-05,-1.84816217627666e-05,-1.66519330021394e-05/

  data alfa2(1,1), alfa2(2,1), alfa2(3,1), alfa2(4,1), alfa2(5,1), &
          alfa2(6,1), alfa2(7,1), alfa2(8,1), alfa2(9,1), alfa2(10,1), &
          alfa2(11,1),alfa2(12,1),alfa2(13,1),alfa2(14,1),alfa2(15,1), &
          alfa2(16,1),alfa2(17,1),alfa2(18,1),alfa2(19,1),alfa2(20,1), &
          alfa2(21,1),alfa2(22,1),alfa2(23,1),alfa2(24,1),alfa2(25,1), &
          alfa2(26,1)     /-3.54211971457744e-04,-1.56161263945159e-04, &
      3.04465503594936e-05, 1.30198655773243e-04, 1.67471106699712e-04, &
      1.70222587683593e-04, 1.56501427608595e-04, 1.36339170977445e-04, &
      1.14886692029825e-04, 9.45869093034688e-05, 7.64498419250898e-05, &
      6.07570334965197e-05, 4.74394299290509e-05, 3.62757512005344e-05, &
      2.69939714979225e-05, 1.93210938247939e-05, 1.30056674793963e-05, &
      7.82620866744497e-06, 3.59257485819352e-06, 1.44040049814252e-07, &
     -2.65396769697939e-06,-4.91346867098486e-06,-6.72739296091248e-06, &
     -8.17269379678658e-06,-9.31304715093561e-06,-1.02011418798016e-05/

  data alfa2(1,2), alfa2(2,2), alfa2(3,2), alfa2(4,2), alfa2(5,2), &
          alfa2(6,2), alfa2(7,2), alfa2(8,2), alfa2(9,2), alfa2(10,2), &
          alfa2(11,2),alfa2(12,2),alfa2(13,2),alfa2(14,2),alfa2(15,2), &
          alfa2(16,2),alfa2(17,2),alfa2(18,2),alfa2(19,2),alfa2(20,2), &
          alfa2(21,2),alfa2(22,2),alfa2(23,2),alfa2(24,2),alfa2(25,2), &
          alfa2(26,2)     / 3.78194199201773e-04, 2.02471952761816e-04, &
     -6.37938506318862e-05,-2.38598230603006e-04,-3.10916256027362e-04, &
     -3.13680115247576e-04,-2.78950273791323e-04,-2.28564082619141e-04, &
     -1.75245280340847e-04,-1.25544063060690e-04,-8.22982872820208e-05, &
     -4.62860730588116e-05,-1.72334302366962e-05, 5.60690482304602e-06, &
      2.31395443148287e-05, 3.62642745856794e-05, 4.58006124490189e-05, &
      5.24595294959114e-05, 5.68396208545815e-05, 5.94349820393104e-05, &
      6.06478527578422e-05, 6.08023907788436e-05, 6.01577894539460e-05, &
      5.89199657344698e-05, 5.72515823777593e-05, 5.52804375585853e-05/

  data beta1(1,1), beta1(2,1), beta1(3,1), beta1(4,1), beta1(5,1), &
          beta1(6,1), beta1(7,1), beta1(8,1), beta1(9,1), beta1(10,1), &
          beta1(11,1),beta1(12,1),beta1(13,1),beta1(14,1),beta1(15,1), &
          beta1(16,1),beta1(17,1),beta1(18,1),beta1(19,1),beta1(20,1), &
          beta1(21,1),beta1(22,1),beta1(23,1),beta1(24,1),beta1(25,1), &
          beta1(26,1)     / 1.79988721413553e-02, 5.59964911064388e-03, &
      2.88501402231133e-03, 1.80096606761054e-03, 1.24753110589199e-03, &
      9.22878876572938e-04, 7.14430421727287e-04, 5.71787281789705e-04, &
      4.69431007606482e-04, 3.93232835462917e-04, 3.34818889318298e-04, &
      2.88952148495752e-04, 2.52211615549573e-04, 2.22280580798883e-04, &
      1.97541838033063e-04, 1.76836855019718e-04, 1.59316899661821e-04, &
      1.44347930197334e-04, 1.31448068119965e-04, 1.20245444949303e-04, &
      1.10449144504599e-04, 1.01828770740567e-04, 9.41998224204238e-05, &
      8.74130545753834e-05, 8.13466262162801e-05, 7.59002269646219e-05/

  data beta1(1,2), beta1(2,2), beta1(3,2), beta1(4,2), beta1(5,2), &
          beta1(6,2), beta1(7,2), beta1(8,2), beta1(9,2), beta1(10,2), &
          beta1(11,2),beta1(12,2),beta1(13,2),beta1(14,2),beta1(15,2), &
          beta1(16,2),beta1(17,2),beta1(18,2),beta1(19,2),beta1(20,2), &
          beta1(21,2),beta1(22,2),beta1(23,2),beta1(24,2),beta1(25,2), &
          beta1(26,2)     /-1.49282953213429e-03,-8.78204709546389e-04, &
     -5.02916549572035e-04,-2.94822138512746e-04,-1.75463996970783e-04, &
     -1.04008550460816e-04,-5.96141953046458e-05,-3.12038929076098e-05, &
     -1.26089735980230e-05,-2.42892608575730e-07, 8.05996165414274e-06, &
      1.36507009262147e-05, 1.73964125472926e-05, 1.98672978842134e-05, &
      2.14463263790823e-05, 2.23954659232457e-05, 2.28967783814713e-05, &
      2.30785389811178e-05, 2.30321976080909e-05, 2.28236073720349e-05, &
      2.25005881105292e-05, 2.20981015361991e-05, 2.16418427448104e-05, &
      2.11507649256221e-05, 2.06388749782171e-05, 2.01165241997082e-05/

  data beta2(1,1), beta2(2,1), beta2(3,1), beta2(4,1), beta2(5,1), &
          beta2(6,1), beta2(7,1), beta2(8,1), beta2(9,1), beta2(10,1), &
          beta2(11,1),beta2(12,1),beta2(13,1),beta2(14,1),beta2(15,1), &
          beta2(16,1),beta2(17,1),beta2(18,1),beta2(19,1),beta2(20,1), &
          beta2(21,1),beta2(22,1),beta2(23,1),beta2(24,1),beta2(25,1), &
          beta2(26,1)     / 5.52213076721293e-04, 4.47932581552385e-04, &
      2.79520653992021e-04, 1.52468156198447e-04, 6.93271105657044e-05, &
      1.76258683069991e-05,-1.35744996343269e-05,-3.17972413350427e-05, &
     -4.18861861696693e-05,-4.69004889379141e-05,-4.87665447413787e-05, &
     -4.87010031186735e-05,-4.74755620890087e-05,-4.55813058138628e-05, &
     -4.33309644511266e-05,-4.09230193157750e-05,-3.84822638603221e-05, &
     -3.60857167535411e-05,-3.37793306123367e-05,-3.15888560772110e-05, &
     -2.95269561750807e-05,-2.75978914828336e-05,-2.58006174666884e-05, &
     -2.41308356761280e-05,-2.25823509518346e-05,-2.11479656768913e-05/

  data beta2(1,2), beta2(2,2), beta2(3,2), beta2(4,2), beta2(5,2), &
          beta2(6,2), beta2(7,2), beta2(8,2), beta2(9,2), beta2(10,2), &
          beta2(11,2),beta2(12,2),beta2(13,2),beta2(14,2),beta2(15,2), &
          beta2(16,2),beta2(17,2),beta2(18,2),beta2(19,2),beta2(20,2), &
          beta2(21,2),beta2(22,2),beta2(23,2),beta2(24,2),beta2(25,2), &
          beta2(26,2)     /-4.74617796559960e-04,-4.77864567147321e-04, &
     -3.20390228067038e-04,-1.61105016119962e-04,-4.25778101285435e-05, &
      3.44571294294968e-05, 7.97092684075675e-05, 1.03138236708272e-04, &
      1.12466775262204e-04, 1.13103642108481e-04, 1.08651634848774e-04, &
      1.01437951597662e-04, 9.29298396593364e-05, 8.40293133016090e-05, &
      7.52727991349134e-05, 6.69632521975731e-05, 5.92564547323195e-05, &
      5.22169308826976e-05, 4.58539485165361e-05, 4.01445513891487e-05, &
      3.50481730031328e-05, 3.05157995034347e-05, 2.64956119950516e-05, &
      2.29363633690998e-05, 1.97893056664022e-05, 1.70091984636413e-05/

  data beta3(1,1), beta3(2,1), beta3(3,1), beta3(4,1), beta3(5,1), &
          beta3(6,1), beta3(7,1), beta3(8,1), beta3(9,1), beta3(10,1), &
          beta3(11,1),beta3(12,1),beta3(13,1),beta3(14,1),beta3(15,1), &
          beta3(16,1),beta3(17,1),beta3(18,1),beta3(19,1),beta3(20,1), &
          beta3(21,1),beta3(22,1),beta3(23,1),beta3(24,1),beta3(25,1), &
          beta3(26,1)     / 7.36465810572578e-04, 8.72790805146194e-04, &
      6.22614862573135e-04, 2.85998154194304e-04, 3.84737672879366e-06, &
     -1.87906003636972e-04,-2.97603646594555e-04,-3.45998126832656e-04, &
     -3.53382470916038e-04,-3.35715635775049e-04,-3.04321124789040e-04, &
     -2.66722723047613e-04,-2.27654214122820e-04,-1.89922611854562e-04, &
     -1.55058918599094e-04,-1.23778240761874e-04,-9.62926147717644e-05, &
     -7.25178327714425e-05,-5.22070028895634e-05,-3.50347750511901e-05, &
     -2.06489761035552e-05,-8.70106096849767e-06, 1.13698686675100e-06, &
      9.16426474122779e-06, 1.56477785428873e-05, 2.08223629482467e-05/

  data gama(1),   gama(2),   gama(3),   gama(4),   gama(5), &
          gama(6),   gama(7),   gama(8),   gama(9),   gama(10), &
          gama(11),  gama(12),  gama(13),  gama(14),  gama(15), &
          gama(16),  gama(17),  gama(18),  gama(19),  gama(20), &
          gama(21),  gama(22),  gama(23),  gama(24),  gama(25), &
          gama(26)        / 6.29960524947437e-01, 2.51984209978975e-01, &
      1.54790300415656e-01, 1.10713062416159e-01, 8.57309395527395e-02, &
      6.97161316958684e-02, 5.86085671893714e-02, 5.04698873536311e-02, &
      4.42600580689155e-02, 3.93720661543510e-02, 3.54283195924455e-02, &
      3.21818857502098e-02, 2.94646240791158e-02, 2.71581677112934e-02, &
      2.51768272973862e-02, 2.34570755306079e-02, 2.19508390134907e-02, &
      2.06210828235646e-02, 1.94388240897881e-02, 1.83810633800683e-02, &
      1.74293213231963e-02, 1.65685837786612e-02, 1.57865285987918e-02, &
      1.50729501494096e-02, 1.44193250839955e-02, 1.38184805735342e-02/
!
  ta = epsilon ( ta )
  tol = max ( ta, 1.0E-15 )
  tb = r1mach(5)
  ju = i1mach(12)

  if ( flgjy /= 1.0E+00 ) then
    jr = i1mach(11)
    elim = 2.303E+00 * tb * ( real ( - ju ) - real ( jr ) )
  else
    elim = 2.303 * ( tb * real ( - ju ) - 3.0E+00 )
  end if

  fn = fnu
  iflw = 0

  do jn = 1, in

    xx = x / fn
    wk(1) = 1.0E+00 - xx * xx
    abw2 = abs ( wk(1) )
    wk(2) = sqrt ( abw2 )
    wk(7) = fn**con2

    if ( abw2 > 0.27750 ) go to 80
!
!  asymptotic expansion
!  cases near x=fn, abs(1.-(x/fn)**2)<=0.2775
!  coefficients of asymptotic expansion by series
!
!  zeta and truncation for a(zeta) and b(zeta) series
!
!  kmax is truncation index for a(zeta) and b(zeta) series=max(2,sa)
!
    if ( abw2 == 0.0E+00 ) then
      sa = 0.0E+00
    else
      sa = tols / log ( abw2 )
    end if

    sb = sa

    do i = 1, 5
      akm = max ( sa, 2.0E+00 )
      kmax(i) = int ( akm )
      sa = sa + sb
    end do

    kb = kmax(5)
    klast = kb - 1
    sa = gama(kb)

    do k = 1, klast
      kb = kb - 1
      sa = sa * wk(1) + gama(kb)
    end do

    z = wk(1) * sa
    az = abs ( z )
    rtz = sqrt ( az )
    wk(3) = con1 * az * rtz
    wk(4) = wk(3) * fn
    wk(5) = rtz * wk(7)
    wk(6) = - wk(5) * wk(5)
    if ( z <= 0.0E+00 ) go to 35
    if ( wk(4) > elim ) go to 75
    wk(6) = - wk(6)

   35   continue

    phi = sqrt ( sqrt ( sa + sa + sa + sa ) )
!
!     b(zeta) for s=0
!
    kb = kmax(5)
    klast = kb - 1
    sb = beta(kb,1)

    do k = 1, klast
      kb = kb - 1
      sb = sb * wk(1) + beta(kb,1)
    end do

    ksp1 = 1
    fn2 = fn * fn
    rfn2 = 1.0E+00 / fn2
    rden = 1.0E+00
    asum = 1.0E+00
    relb = tol * abs ( sb )
    bsum = sb

    do ks = 1, 4

      ksp1 = ksp1 + 1
      rden = rden * rfn2
!
!     a(zeta) and b(zeta) for s=1,2,3,4
!
      kstemp = 5 - ks
      kb = kmax(kstemp)
      klast = kb - 1
      sa = alfa(kb,ks)
      sb = beta(kb,ksp1)

      do k = 1, klast
        kb = kb - 1
        sa = sa * wk(1) + alfa(kb,ks)
        sb = sb * wk(1) + beta(kb,ksp1)
      end do

      ta = sa * rden
      tb = sb * rden
      asum = asum + ta
      bsum = bsum + tb

      if ( abs ( ta ) <= tol .and. abs ( tb ) <= relb ) then
        exit
      end if

    end do

!  70   continue

    bsum = bsum / ( fn * wk(7) )

    go to 160

   75   continue
    iflw = 1
    return

   80   continue

    upol(1) = 1.0E+00
    tau = 1.0E+00 / wk(2)
    t2 = 1.0E+00 / wk(1)
    if ( wk(1) >= 0.0E+00 ) go to 90
!
!  cases for (x/fn)>sqrt(1.2775)
!
    wk(3) = abs(wk(2)-atan(wk(2)))
    wk(4) = wk(3)*fn
    rcz = -con1/wk(4)
    z32 = 1.5*wk(3)
    rtz = z32**con2
    wk(5) = rtz*wk(7)
    wk(6) = -wk(5)*wk(5)
    go to 100
   90   continue
!
!     cases for (x/fn)<sqrt(0.7225)
!
    wk(3) = abs( log((1.0+wk(2))/xx)-wk(2))
    wk(4) = wk(3)*fn
    rcz = con1/wk(4)
    if ( wk(4)>elim) go to 75
    z32 = 1.5*wk(3)
    rtz = z32**con2
    wk(7) = fn**con2
    wk(5) = rtz*wk(7)
    wk(6) = wk(5)*wk(5)
  100   continue
    phi = sqrt((rtz+rtz)*tau)
    tb = 1.0E+00
    asum = 1.0E+00
    tfn = tau/fn
    upol(2) = (c(1)*t2+c(2))*tfn
    crz32 = con548*rcz
    bsum = upol(2) + crz32
    relb = tol*abs(bsum)
    ap = tfn
    ks = 0
    kp1 = 2
    rzden = rcz
    l = 2

    do lr=2,8,2
!
!     compute two u polynomials for next a(zeta) and b(zeta)
!
      lrp1 = lr + 1
      do k=lr,lrp1
        ks = ks + 1
        kp1 = kp1 + 1
        l = l + 1
        s1 = c(l)
        do j=2,kp1
          l = l + 1
          s1 = s1*t2 + c(l)
        end do
        ap = ap*tfn
        upol(kp1) = ap*s1
        cr(ks) = br(ks)*rzden
        rzden = rzden*rcz
        dr(ks) = ar(ks)*rzden
      end do

      suma = upol(lrp1)
      sumb = upol(lr+2) + upol(lrp1)*crz32
      ju = lrp1

      do jr = 1, lr
        ju = ju - 1
        suma = suma + cr(jr)*upol(ju)
        sumb = sumb + dr(jr)*upol(ju)
      end do

      tb = -tb
      if (wk(1)>0.0) tb = abs(tb)
      asum = asum + suma*tb
      bsum = bsum + sumb*tb

      if ( abs(suma)<=tol .and. abs(sumb)<=relb ) then
        exit
      end if

    end do

! 150 continue

    tb = wk(5)
    if ( wk(1) > 0.0E+00 ) tb = -tb
    bsum = bsum/tb

  160   continue
    call funjy(wk(6), wk(5), wk(4), fi, dfi)
    y(jn) = flgjy*phi*(fi*asum+dfi*bsum)/wk(7)
    fn = fn - flgjy

  end do

  return
end
subroutine bakslv ( nr, n, a, x, b )
!
!*******************************************************************************
!
!! BAKSLV solves A*x=b where A is an upper triangular matrix.
!
!
!  Discussion:
!
!    BAKSLV solves the linear equations TRANSPOSE(A)*x=b, where A is a
!    lower triangular matrix.  This routine is required by the UNCMIN
!    minimization program.
!
!    If B is no longer required by calling routine, then vectors B and
!    X may share the same storage, and the output value of X will
!    overwrite the input value of B.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, integer NR, the leading dimension of A.
!
!    Input, integer N, the number of rows and columns in A.
!
!    Input, real A(NR,N), the N by N matrix, containing the lower
!    triangular matrix.  A is not altered by this routine.
!
!    Output, real X(N), the solution vector.
!
!    Input, real B(N), the right hand side vector.
!
  integer n
  integer nr
!
  real a(nr,n)
  real b(n)
  integer i
  integer ip1
  integer j
  real x(n)
!
!  Solve (l-transpose)x=b. (back solve)
!
  i = n
  x(i) = b(i) / a(i,i)

  if ( n == 1 ) then
    return
  end if

  do

    ip1 = i
    i = i - 1

    x(i) = ( b(i) - dot_product ( x(ip1:n), a(ip1:n,i) ) ) / a(i,i)

    if ( i == 1 ) then
      exit
    end if

  end do

  return
end
function besi0 ( x )
!
!*******************************************************************************
!
!! BESI0 computes the hyperbolic Bessel function of the first kind of order zero.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, real X, the argument of the Bessel function.
!
!    Output, real BESI0, the value of the Bessel function at X.
!
  real besi0
  real besi0e
  real bi0cs(12)
  real csevl
  integer inits
  integer, save :: nti0 = 0
  real r1mach
  real x
  real, save :: xmax = 0.0E+00
  real, save :: xsml = 0.0E+00
  real y
!
  data bi0cs( 1) /   -0.07660547252839144951E+00 /
  data bi0cs( 2) /    1.927337953993808270E+00 /
  data bi0cs( 3) /    0.2282644586920301339E+00 /
  data bi0cs( 4) /    0.01304891466707290428E+00 /
  data bi0cs( 5) /    0.00043442709008164874E+00 /
  data bi0cs( 6) /    0.00000942265768600193E+00 /
  data bi0cs( 7) /    0.00000014340062895106E+00 /
  data bi0cs( 8) /    0.00000000161384906966E+00 /
  data bi0cs( 9) /    0.00000000001396650044E+00 /
  data bi0cs(10) /    0.00000000000009579451E+00 /
  data bi0cs(11) /    0.00000000000000053339E+00 /
  data bi0cs(12) /    0.00000000000000000245E+00 /
!
  if ( nti0 == 0 ) then
    nti0 = inits ( bi0cs, 12, 0.1E+00 * r1mach(3) )
    xsml = 2.0E+00 * sqrt ( epsilon ( xsml ) )
    xmax = log ( huge ( xmax ) )
  end if

  y = abs ( x )

  if ( y <= 3.0E+00 ) then

    if ( y > xsml ) then
      besi0 = 2.75E+00 + csevl ( y*y/4.5E+00-1.0E+00, bi0cs, nti0 )
    else
      besi0 = 1.0E+00
    end if

    return

  end if

  if ( y > xmax ) then
    write ( *, * ) ' '
    write ( *, * ) 'BESI0 - Fatal error!'
    write ( *, * ) '  |X| is so big that BESI0 will overflow.'
    stop
  end if

  besi0 = exp(y) * besi0e(x)

  return
end
function besi0e ( x )
!
!*******************************************************************************
!
!! BESI0E computes the exponentially scaled hyperbolic Bessel function of the first kind of order ze
!
!
!  Discussion:
!
!    BESI0E calculates the exponentially scaled modified hyperbolic
!    Bessel function of the first kind of order zero for real argument X.
!
!      besi0e(x) = exp ( - abs ( x ) ) * i0 ( x ).
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, real X, the argument of the Bessel function.
!
!    Output, real BESI0E, the value of the Bessel function at X.
!
  real ai02cs(22)
  real ai0cs(21)
  real besi0e
  real bi0cs(12)
  real csevl
  integer inits
  integer ntai0
  integer ntai02
  integer nti0
  real r1mach
  real x
  real, save :: xsml = 0.0E+00
  real y
!
  data bi0cs( 1) /   -.07660547252839144951 /
  data bi0cs( 2) /   1.927337953993808270 /
  data bi0cs( 3) /    .2282644586920301339 /
  data bi0cs( 4) /    .01304891466707290428 /
  data bi0cs( 5) /    .00043442709008164874 /
  data bi0cs( 6) /    .00000942265768600193 /
  data bi0cs( 7) /    .00000014340062895106 /
  data bi0cs( 8) /    .00000000161384906966 /
  data bi0cs( 9) /    .00000000001396650044 /
  data bi0cs(10) /    .00000000000009579451 /
  data bi0cs(11) /    .00000000000000053339 /
  data bi0cs(12) /    .00000000000000000245 /
  data ai0cs( 1) /    .07575994494023796 /
  data ai0cs( 2) /    .00759138081082334 /
  data ai0cs( 3) /    .00041531313389237 /
  data ai0cs( 4) /    .00001070076463439 /
  data ai0cs( 5) /   -.00000790117997921 /
  data ai0cs( 6) /   -.00000078261435014 /
  data ai0cs( 7) /    .00000027838499429 /
  data ai0cs( 8) /    .00000000825247260 /
  data ai0cs( 9) /   -.00000001204463945 /
  data ai0cs(10) /    .00000000155964859 /
  data ai0cs(11) /    .00000000022925563 /
  data ai0cs(12) /   -.00000000011916228 /
  data ai0cs(13) /    .00000000001757854 /
  data ai0cs(14) /    .00000000000112822 /
  data ai0cs(15) /   -.00000000000114684 /
  data ai0cs(16) /    .00000000000027155 /
  data ai0cs(17) /   -.00000000000002415 /
  data ai0cs(18) /   -.00000000000000608 /
  data ai0cs(19) /    .00000000000000314 /
  data ai0cs(20) /   -.00000000000000071 /
  data ai0cs(21) /    .00000000000000007 /
  data ai02cs( 1) /    .05449041101410882 /
  data ai02cs( 2) /    .00336911647825569 /
  data ai02cs( 3) /    .00006889758346918 /
  data ai02cs( 4) /    .00000289137052082 /
  data ai02cs( 5) /    .00000020489185893 /
  data ai02cs( 6) /    .00000002266668991 /
  data ai02cs( 7) /    .00000000339623203 /
  data ai02cs( 8) /    .00000000049406022 /
  data ai02cs( 9) /    .00000000001188914 /
  data ai02cs(10) /   -.00000000003149915 /
  data ai02cs(11) /   -.00000000001321580 /
  data ai02cs(12) /   -.00000000000179419 /
  data ai02cs(13) /    .00000000000071801 /
  data ai02cs(14) /    .00000000000038529 /
  data ai02cs(15) /    .00000000000001539 /
  data ai02cs(16) /   -.00000000000004151 /
  data ai02cs(17) /   -.00000000000000954 /
  data ai02cs(18) /    .00000000000000382 /
  data ai02cs(19) /    .00000000000000176 /
  data ai02cs(20) /   -.00000000000000034 /
  data ai02cs(21) /   -.00000000000000027 /
  data ai02cs(22) /    .00000000000000003 /
  data nti0, ntai0, ntai02 / 3*0 /
!
  if ( nti0 == 0 ) then
    nti0 = inits ( bi0cs, 12, 0.1E+00*r1mach(3) )
    ntai0 = inits ( ai0cs, 21, 0.1E+00*r1mach(3) )
    ntai02 = inits ( ai02cs, 22, 0.1E+00*r1mach(3) )
    xsml = 2.0E+00 * sqrt ( epsilon ( xsml ) )
  end if

  y = abs(x)

       if ( y <= xsml ) then

    besi0e = 1.0E+00

  else if ( y <= 3.0E+00 ) then

    besi0e = exp(-y) * ( 2.75 + csevl (y*y/4.5-1.0, bi0cs, nti0) )

  else if ( y <= 8.0E+00 ) then

    besi0e = ( 0.375E+00 + csevl ((48.0/y-11.0)/5.0, ai0cs, ntai0) ) / sqrt(y)

  else if ( y > 8.0E+00 ) then

    besi0e = ( 0.375E+00 + csevl (16.0/y-1.0, ai02cs, ntai02)) / sqrt(y)

  end if

  return
end
subroutine besj ( x, alpha, n, y, nz )
!
!*******************************************************************************
!
!! BESJ computes an N member sequence of J Bessel functions.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    D. E. Amos, S. L. Daniel, M. K. Weston,
!    SAND-75-0147
!    CDC 6600 subroutines IBESS and JBESS for Bessel functions
!      I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0
!    ACM Transactions on Mathematical Software,
!    Volume 3, pages 76-92, 1977.
!
!    F. W. J. Olver,
!    Tables of Bessel Functions of Moderate or Large Orders,
!    NPL Mathematical Tables, Vol. 6,
!    Her Majesty's Stationery Office, London, 1962.
!
!  Discussion:
!
!    BESJ computes an N member sequence of J Bessel functions
!
!      J(ALPHA+K-1) (X)
!
!    for K=1,..,N for non-negative order ALPHA and X.
!
!    A combination of the power series, the asymptotic expansion for X
!    to infinity and the uniform asymptotic expansion for NU to infinity
!    are applied over subdivisions of the (NU,X) plane.  For values of
!    (NU,X) not covered by one of these formulae, the order is
!    incremented or decremented by integer values into a region where
!    one of the formulas apply.
!
!    Backward recursion is applied to reduce orders by integer values
!    except where the entire sequence lies in the oscillatory region.
!    In this case forward recursion is stable and values from the
!    asymptotic expansion for X to infinity start the recursion when it
!    is efficient to do so.
!
!    Leading terms of the series and uniform expansion are tested for
!    underflow.  If a sequence is requested and the last member would
!    underflow, the result is set to zero and the next lower order
!    tried, etc., until a member comes on scale or all members are set
!    to zero.
!
!    Overflow cannot occur.
!
!  Parameters:
!
!    Input, real X, X .GE. 0.0, the argument of the Bessel function.
!
!    Input, real ALPHA, order of first member of the sequence.
!    ALPHA must be at least 0.0.
!
!    Input, integer N, number of members in the sequence,
!    N must be at least 1.
!
!    Output, real Y(N), a vector whose first N components contain
!    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
!
!    Output, integer NZ, the number of components of Y set to zero
!    due to underflow,
!
!    NZ=0, normal return, computation completed
!
!    NZ .NE. 0, Y(N-NZ+1) through Y(N) were set to 0.
!
  external jairy
!
  real ak
  real akm
  real alngam
  real alpha
  real ans
  real ap
  real arg
  real coef
  real dalpha
  real dfn
  real dtm
  real earg
  real elim1
  real etx
  real fidal
  real flgjy
  real fn
  real fnf
  real fni
  real fnp1
  real fnu
  real fnulim(2)
  real gln
  integer i
  integer i1
  integer i1mach
  integer i2
  integer ialp
  integer idalp
  integer iflw
  integer in
  integer, parameter :: inlim = 150
  integer is
  integer k
  integer kk
  integer km
  integer kt
  integer n
  integer nn
  integer ns
  integer nz
  real, parameter :: pdf = 0.785398163397448E+00
  real pidt
  real pp(4)
  real r1mach
  real rden
  real relb
  real rttp
  real, parameter :: rtwo = 1.34839972492648E+00
  real rtx
  real rzden
  real s
  real sa
  real sb
  real sxo2
  real s1
  real s2
  real t
  real ta
  real tau
  real tb
  real temp(3)
  real tfn
  real tm
  real tol
  real tolln
  real trx
  real tx
  real t1
  real t2
  real wk(7)
  real x
  real xo2
  real xo2l
  real y(1)
!
  data rttp,pidt / 7.97884560802865e-01, 1.57079632679490/
  data  pp(1),  pp(2),  pp(3),  pp(4) / 8.72909153935547, &
    2.65693932265030e-01, 1.24578576865586e-01, 7.70133747430388e-04/
  data fnulim(1), fnulim(2) /      100.0,     60.0E+00     /
!
  nz = 0
  kt = 1
!     i1mach(14) replaces i1mach(11) in a double precision code
!     i1mach(15) replaces i1mach(12) in a double precision code
  ta = epsilon ( ta )
  tol = max (ta,1.0e-15)
  i1 = i1mach(11) + 1
  i2 = i1mach(12)
  tb = r1mach(5)
  elim1 = 2.303*(real(-i2)*tb-3.0)
!
!     tolln = -ln(tol)
!
  tolln = 2.303*tb*real(i1)
  tolln = min (tolln,34.5388)
  if (n-1) 720, 10, 20
   10 kt = 2
   20 nn = n
  if (x) 730, 30, 80
   30 if (alpha) 710, 40, 50
   40 y(1) = 1.0E+00
  if (n==1) return
  i1 = 2
  go to 60
   50 i1 = 1
   60 continue

  y(i1:n) = 0.0E+00

  return
   80 continue
  if (alpha<0.0) go to 710

  ialp = int(alpha)
  fni = real(ialp+n-1)
  fnf = alpha - real(ialp)
  dfn = fni + fnf
  fnu = dfn
  xo2 = x*0.5
  sxo2 = xo2*xo2
!
!     decision tree for region where series, asymptotic expansion for x
!     to infinity and asymptotic expansion for nu to infinity are
!     applied.
!
  if (sxo2<=(fnu+1.0)) go to 90
  ta = max (20.0,fnu)
  if (x>ta) go to 120
  if (x>12.0) go to 110
  xo2l = log(xo2)
  ns = int(sxo2-fnu) + 1
  go to 100
   90 fn = fnu
  fnp1 = fn + 1.0E+00
  xo2l = log(xo2)
  is = kt
  if (x<=0.50) go to 330
  ns = 0
  100 fni = fni + real(ns)
  dfn = fni + fnf
  fn = dfn
  fnp1 = fn + 1.0E+00
  is = kt
  if (n-1+ns>0) is = 3
  go to 330
  110 ans = max (36.0-fnu,0.0)
  ns = int(ans)
  fni = fni + real(ns)
  dfn = fni + fnf
  fn = dfn
  is = kt
  if (n-1+ns>0) is = 3
  go to 130
  120 continue
  rtx = sqrt(x)
  tau = rtwo*rtx
  ta = tau + fnulim(kt)
  if (fnu<=ta) go to 480
  fn = fnu
  is = kt
!
!  Uniform asymptotic expansion for nu to infinity
!
  130 continue
  i1 = abs(3-is)
  i1 = max(i1,1)
  flgjy = 1.0E+00
  call asyjy(jairy,x,fn,flgjy,i1,temp(is),wk,iflw)
  if ( iflw/=0) go to 380
  go to (320, 450, 620), is
  310 temp(1) = temp(3)
  kt = 1
  320 is = 2
  fni = fni - 1.0E+00
  dfn = fni + fnf
  fn = dfn
  if ( i1==2) go to 450
  go to 130
!
!  Series for (x/2)**2<=nu+1
!
  330 continue
  gln = alngam(fnp1)
  arg = fn*xo2l - gln
  if (arg<(-elim1)) go to 400
  earg = exp(arg)
  340 continue
  s = 1.0E+00
  if (x<tol) go to 360
  ak = 3.0E+00
  t2 = 1.0E+00
  t = 1.0E+00
  s1 = fn

  do k=1, 17
    s2 = t2 + s1
    t = -t*sxo2/s2
    s = s + t
    if (abs(t)<tol) go to 360
    t2 = t2 + ak
    ak = ak + 2.0E+00
    s1 = s1 + fn
  end do

  360 continue
  temp(is) = s*earg
  go to (370, 450, 610), is
  370 earg = earg*fn/xo2
  fni = fni - 1.0E+00
  dfn = fni + fnf
  fn = dfn
  is = 2
  go to 340
!
!  Set underflow value and update parameters
!
  380 y(nn) = 0.0E+00
  nn = nn - 1
  fni = fni - 1.0E+00
  dfn = fni + fnf
  fn = dfn
  if (nn-1) 440, 390, 130
  390 kt = 2
  is = 2
  go to 130
  400 y(nn) = 0.0E+00
  nn = nn - 1
  fnp1 = fn
  fni = fni - 1.0E+00
  dfn = fni + fnf
  fn = dfn
  if (nn-1) 440, 410, 420
  410 continue
  kt = 2
  is = 2
  420 if (sxo2<=fnp1) go to 430
  go to 130
  430 arg = arg - xo2l + log(fnp1)
  if (arg<(-elim1)) go to 400
  go to 330
  440 nz = n - nn
  return
!
!  Backward recursion section
!
  450 continue
  nz = n - nn
  if (kt==2) go to 470
!     backward recur from index alpha+nn-1 to alpha
  y(nn) = temp(1)
  y(nn-1) = temp(2)
  if (nn==2) return
  trx = 2.0/x
  dtm = fni
  tm = (dtm+fnf)*trx
  k = nn + 1
  do i=3,nn
    k = k - 1
    y(k-2) = tm*y(k-1) - y(k)
    dtm = dtm - 1.0E+00
    tm = (dtm+fnf)*trx
  end do
  return
  470 y(1) = temp(2)
  return
!
!     asymptotic expansion for x to infinity with forward recursion in
!     oscillatory region x>max(20, nu), provided the last member
!     of the sequence is also in the region.
!
  480 continue
  in = int(alpha-tau+2.0)
  if (in<=0) go to 490
  idalp = ialp - in - 1
  kt = 1
  go to 500
  490 continue
  idalp = ialp
  in = 0
  500 is = kt
  fidal = real(idalp)
  dalpha = fidal + fnf
  arg = x - pidt*dalpha - pdf
  sa = sin(arg)
  sb = cos(arg)
  coef = rttp/rtx
  etx = 8.0*x
  510 continue
  dtm = fidal + fidal
  dtm = dtm*dtm
  tm = 0.0E+00
  if (fidal==0.0E+00 .and. abs(fnf)<tol) go to 520
  tm = 4.0*fnf*(fidal+fidal+fnf)
  520 continue
  trx = dtm - 1.0E+00
  t2 = (trx+tm)/etx
  s2 = t2
  relb = tol*abs(t2)
  t1 = etx
  s1 = 1.0E+00
  fn = 1.0E+00
  ak = 8.0E+00

  do k=1,13
    t1 = t1 + etx
    fn = fn + ak
    trx = dtm - fn
    ap = trx + tm
    t2 = -t2*ap/t1
    s1 = s1 + t2
    t1 = t1 + etx
    ak = ak + 8.0E+00
    fn = fn + ak
    trx = dtm - fn
    ap = trx + tm
    t2 = t2*ap/t1
    s2 = s2 + t2
    if (abs(t2)<=relb) go to 540
    ak = ak + 8.0E+00
  end do

  540 temp(is) = coef*(s1*sb-s2*sa)
  if ( is==2) go to 560
  550 fidal = fidal + 1.0E+00
  dalpha = fidal + fnf
  is = 2
  tb = sa
  sa = -sb
  sb = tb
  go to 510
!
!  Forward recursion section
!
  560 if (kt==2) go to 470
  s1 = temp(1)
  s2 = temp(2)
  tx = 2.0/x
  tm = dalpha*tx
  if (in==0) go to 580
!
!   Forward recur to index alpha
!
  do i = 1, in
    s = s2
    s2 = tm*s2 - s1
    tm = tm + tx
    s1 = s
  end do

  if (nn==1) go to 600
  s = s2
  s2 = tm*s2 - s1
  tm = tm + tx
  s1 = s
  580 continue
!
!  Forward recur from index alpha to alpha+n-1
!
  y(1) = s1
  y(2) = s2
  if (nn==2) return
  do i=3,nn
    y(i) = tm*y(i-1) - y(i-2)
    tm = tm + tx
  end do
  return
  600 y(1) = s2
  return
!
!     backward recursion with normalization by
!     asymptotic expansion for nu to infinity or power series.
!
  610 continue
!     computation of last order for series normalization
  akm = max (3.0-fn,0.0)
  km = int(akm)
  tfn = fn + real(km)
  ta = (gln+tfn-0.9189385332-0.0833333333/tfn)/(tfn+0.5)
  ta = xo2l - ta
  tb = -(1.0-1.5/tfn)/tfn
  akm = tolln/(-ta+sqrt(ta*ta-tolln*tb)) + 1.5
  in = km + int(akm)
  go to 660
  620 continue
!
!     computation of last order for asymptotic expansion normalization
!
  gln = wk(3) + wk(2)
  if (wk(6)>30.0) go to 640
  rden = (pp(4)*wk(6)+pp(3))*wk(6) + 1.0E+00
  rzden = pp(1) + pp(2)*wk(6)
  ta = rzden/rden
  if (wk(1)<0.10) go to 630
  tb = gln/wk(5)
  go to 650
  630 continue

  tb=(1.259921049+(0.1679894730+0.0887944358*wk(1))*wk(1))/wk(7)
  go to 650
  640 continue
  ta = 0.5*tolln/wk(4)
  ta=((0.0493827160*ta-0.1111111111)*ta+0.6666666667)*ta*wk(6)
  if (wk(1)<0.10) go to 630
  tb = gln/wk(5)
  650 in = int(ta/tb+1.5)
  if (in>inlim) go to 310
  660 continue
  dtm = fni + real(in)
  trx = 2.0/x
  tm = (dtm+fnf)*trx
  ta = 0.0E+00
  tb = tol
  kk = 1
  670 continue
!
!     backward recur unindexed
!
  do i=1,in
    s = tb
    tb = tm*tb - ta
    ta = s
    dtm = dtm - 1.0E+00
    tm = (dtm+fnf)*trx
  end do
!
!     normalization
  if (kk/=1) go to 690
  ta = (ta/tb)*temp(3)
  tb = temp(3)
  kk = 2
  in = ns
  if (ns/=0) go to 670
  690 y(nn) = tb
  nz = n - nn
  if (nn==1) return
  k = nn - 1
  y(k) = tm*tb - ta
  if (nn==2) return
  dtm = dtm - 1.0E+00
  tm = (dtm+fnf)*trx
  km = k - 1
!
!  backward recur indexed
!
  do i=1,km
    y(k-1) = tm*y(k) - y(k+1)
    dtm = dtm - 1.0E+00
    tm = (dtm+fnf)*trx
    k = k - 1
  end do

  return

  710 continue
  call xerror( 'besj - order, alpha, less than zero.', 36, 2, 1)
  return
  720 continue
  call xerror( 'besj - n less than one.', 23, 2, 1)
  return
  730 continue
  call xerror( 'besj - x less than zero.', 24, 2, 1)
  return
end
subroutine bp ( n, b, x )
!
!*******************************************************************************
!
!! BP evaluates the N+1 Bernstein basis functions of degree N on [0,1].
!
!
!  Definition:
!
!    The I-th Bernstein basis polynomial of degree N is defined as:
!
!      B(N,I,X)= N!/(I!*(N-I)!) * (1-X)**(N-I) * X**I
!
!    although this is not how the values are computed.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, integer N, should be 0 or greater.
!
!    Input, real X, the point where the functions should be evaluated.
!
!    Output, real B(0:N), the values of the Bernstein polynomials
!    at the point X.
!
  integer n
!
  real b(0:n)
  integer i
  integer j
  real x
!
  if ( n == 0 ) then

    b(0) = 1.0E+00

  else if ( n > 0 ) then

    do i = 1, n

       if ( i == 1 ) then
         b(1) = x
       else
         b(i) = x * b(i-1)
       end if

       do j = i-1, 1, -1
         b(j) = x * b(j-1) + ( 1.0E+00 - x ) * b(j)
       end do

       if ( i == 1 ) then
         b(0) = 1.0E+00 - x
       else
         b(0) = ( 1.0E+00 - x ) * b(0)
       end if

    end do

  end if

  return
end
subroutine cfftb ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTB computes the backward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    CFFTB computes a complex periodic sequence from its Fourier coefficients.
!
!    A call of CFFTF followed by a call of CFFTB will multiply the
!    sequence by N.  In other words, the transforms are not normalized.
!
!    The array WSAVE must be initialized by CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N )
!        C_in(K) * exp ( sqrt ( - 1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N.
!
  integer n
!
  complex c(n)
  real wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftb1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftb1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTB1 is a lower-level routine used by CFFTB.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input/output, complex C(N).
!    On input, C contains the sequence of Fourier coefficients.
!    On output, C contains the sequence of data values that correspond
!    to the input coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  complex c(n)
  complex ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nac
  integer nf
  real wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passb4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passb4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passb2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passb2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passb3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passb3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passb5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passb5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passb ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passb ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cfftb_2d ( ldf, n, f, wsave )
!
!*******************************************************************************
!
!! CFFTB_2D computes a backward two dimensional complex fast Fourier transform.
!
!
!  Discussion:
!
!    The routine computes the backward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to CFFTB_2D followed by a call
!    to CFFTF_2D will return the original data multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) * (-1.)**(I+J),  I,J =0,...,N-1.
!
!    Before calling CFFTF_2D or CFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling CFFTI.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Modified:
!
!    12 March 2001
!
!  Parameters:
!
!    Input, integer LDF, the leading dimension of the matrix.
!
!    Input, integer N, the number of rows and columns in the matrix.
!
!    Input/output, complex F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real WSAVE(4*N+15), a work array whose values depend on N,
!    and which must be initialized by calling CFFTI.
!
  integer ldf
  integer n
!
  complex f(ldf,n)
  integer i
  real wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call cfftb ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call cfftb ( n, f(1,i), wsave )
  end do

  return
end
subroutine cfftf ( n, c, wsave )
!
!*******************************************************************************
!
!! CFFTF computes the forward complex discrete Fourier transform.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    CFFTF computes the Fourier coefficients of a complex periodic sequence.
!
!    The transform is not normalized.  To obtain a normalized transform,
!    the output must be divided by N.  Otherwise a call of CFFTF
!    followed by a call of CFFTB will multiply the sequence by N.
!
!    The array WSAVE must be initialized by calling CFFTI.
!
!    The transform is defined by:
!
!      C_out(J) = sum ( 1 <= K <= N )
!        C_in(K) * exp ( - sqrt ( -1 ) * ( J - 1 ) * ( K - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is more efficient when N is the product of small primes.
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, real WSAVE(4*N+15).  The array must be initialized by calling
!    CFFTI.  A different WSAVE array must be used for each different
!    value of N.
!
  integer n
!
  complex c(n)
  real wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cfftf1 ( n, c, wsave(1), wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! CFFTF1 is a lower level routine used by CFFTF.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input/output, complex C(N).
!    On input, the data sequence to be transformed.
!    On output, the Fourier coefficients.
!
!    Input, complex CH(N).
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  complex c(n)
  complex ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nac
  integer nf
  real wa(2*n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = 2 * ido * l1

    if ( ip == 4 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido

      if ( na == 0 ) then
        call passf4 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call passf4 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call passf2 ( 2*ido, l1, c, ch, wa(iw) )
      else
        call passf2 ( 2*ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + 2 * ido

      if ( na == 0 ) then
        call passf3 ( 2*ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call passf3 ( 2*ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + 2 * ido
      ix3 = ix2 + 2 * ido
      ix4 = ix3 + 2 * ido

      if ( na == 0 ) then
        call passf5 ( 2*ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call passf5 ( 2*ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call passf ( nac, 2*ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call passf ( nac, 2*ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( nac /= 0 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * 2 * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine cfftf_2d ( ldf, n, f, wsave )
!
!*******************************************************************************
!
!! CFFTF_2D computes a two dimensional complex fast Fourier transform.
!
!
!  Discussion:
!
!    The routine computes the forward two dimensional fast Fourier transform,
!    of a complex N by N matrix of data.
!
!    The output is unscaled, that is, a call to CFFTF_2D,
!    followed by a call to CFFTB_2D will return the original data
!    multiplied by N*N.
!
!    For some applications it is desirable to have the transform scaled so
!    the center of the N by N frequency square corresponds to zero
!    frequency.  The user can do this replacing the original input data
!    F(I,J) by F(I,J) *(-1.)**(I+J),  I,J =0,...,N-1.
!
!    Before calling CFFTF_2D or CFFTB_2D, it is necessary to initialize
!    the array WSAVE by calling CFFTI.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Modified:
!
!    12 March 2001
!
!  Parameters:
!
!    Input, integer LDF, the leading dimension of the matrix.
!
!    Input, integer N, the number of rows and columns in the matrix.
!
!    Input/output, complex F(LDF,N),
!    On input, an N by N array of complex values to be transformed.
!    On output, the transformed values.
!
!    Input, real WSAVE(4*N+15), a work array whose values depend on N,
!    and which must be initialized by calling CFFTI.
!
  integer ldf
  integer n
!
  complex f(ldf,n)
  integer i
  real wsave(4*n+15)
!
!  Row transforms:
!
  f(1:n,1:n) = transpose ( f(1:n,1:n) )

  do i = 1, n
    call cfftf ( n, f(1,i), wsave )
  end do

  f(1:n,1:n) = transpose ( f(1:n,1:n) )
!
!  Column transforms:
!
  do i = 1, n
    call cfftf ( n, f(1,i), wsave )
  end do

  return
end
subroutine cffti ( n, wsave )
!
!*******************************************************************************
!
!! CFFTI initializes WSAVE, used in CFFTF and CFFTB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Output, real WSAVE(4*N+15), contains data, dependent on the value
!    of N, which is necessary for the CFFTF or CFFTB routines.
!
  integer n
!
  real wsave(4*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call cffti1 ( n, wsave(2*n+1), wsave(4*n+1) )

  return
end
subroutine cffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! CFFTI1 is a lower level routine used by CFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input, real WA(2*N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  real arg
  real argh
  real argld
  real fi
  integer i
  integer i1
  integer ib
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer j
  integer k1
  integer l1
  integer l2
  integer ld
  integer nf
  real pimach
  real wa(2*n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0E+00 * pimach() / real ( n )
  i = 2
  l1 = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      i1 = i
      wa(i-1) = 1.0E+00
      wa(i) = 0.0E+00
      ld = ld + l1
      fi = 0.0E+00
      argld = real ( ld ) * argh

      do ii = 4, 2*ido+2, 2
        i = i + 2
        fi = fi + 1.0E+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      if ( ip > 5 ) then
        wa(i1-1) = wa(i-1)
        wa(i1) = wa(i)
      end if

    end do

    l1 = l2

  end do

  return
end
subroutine chfdv ( x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr )
!
!*******************************************************************************
!
!! CHFDV evaluates a cubic polynomial and its derivative given in Hermite form.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Discussion:
!
!    CHFDV evaluates a cubic polynomial given in Hermite form and its
!    first derivative at an array of points.  While designed for
!    use by PCHFD, it may be useful directly as an evaluator for
!    a piecewise cubic Hermite function in applications, such as
!    graphing, where the interval is known in advance.
!
!    If only function values are required, use CHFEV instead.
!
!  Parameters:
!
!     x1,x2 -- (input) endpoints of interval of definition of cubic.
!           (error return if  x1==x2 .)
!
!     f1,f2 -- (input) values of function at x1 and x2, respectively.
!
!     d1,d2 -- (input) values of derivative at x1 and x2, respectively.
!
!     ne -- (input) number of evaluation points.  (error return if
!           ne<1 .)
!
!     xe -- (input) real array of points at which the functions are to
!           be evaluated.  if any of the xe are outside the interval
!           [x1,x2], a warning error is returned in next.
!
!     fe -- (output) real array of values of the cubic function defined
!           by  x1,x2, f1,f2, d1,d2  at the points  xe.
!
!     de -- (output) real array of values of the first derivative of
!           the same function at the points  xe.
!
!     next -- (output) integer array indicating number of extrapolation
!           points:
!            next(1) = number of evaluation points to left of interval.
!            next(2) = number of evaluation points to right of interval.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           "recoverable" errors:
!              ierr = -1  if ne<1 .
!              ierr = -2  if x1==x2 .
!                (output arrays have not been changed in either case.)
!
  integer ne
!
  real c2
  real c2t2
  real c3
  real c3t3
  real d1
  real d2
  real de(ne)
  real del1
  real del2
  real delta
  real f1
  real f2
  real fe(ne)
  real h
  integer i
  integer ierr
  integer next(2)
  real x
  real x1
  real x2
  real xe(ne)
  real xma
  real xmi
!
!  validity-check arguments.
!
  if ( ne < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CHFDV - Fatal error!'
    write ( *, * ) '  The number of evaluation points was less than 1.'
    stop
  end if

  h = x2 - x1
  if (h == 0.0)  go to 5002
!
!  initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min (0.0, h)
  xma = max (0.0, h)
!
!  compute cubic coefficients (expanded about x1).
!
  delta = (f2 - f1)/h
  del1 = (d1 - delta)/h
  del2 = (d2 - delta)/h
!                                           (delta is no longer needed.)
  c2 = -(del1+del1 + del2)
  c2t2 = c2 + c2
  c3 = (del1 + del2)/h
!                               (h, del1 and del2 are no longer needed.)
  c3t3 = c3+c3+c3
!
!  evaluation loop.
!
  do i = 1, ne
     x = xe(i) - x1
     fe(i) = f1 + x*(d1 + x*(c2 + x*c3))
     de(i) = d1 + x*(c2t2 + x*c3t3)
!          count extrapolation points.
     if ( x<xmi )  next(1) = next(1) + 1
     if ( x>xma )  next(2) = next(2) + 1
!        (note redundancy--if either condition is true, other is false.)
  end do

  return
!
!  error returns.
!
 5002 continue

  ierr = -2
  call xerror ('chfdv -- interval endpoints equal', 33, ierr, 1)
  return
end
subroutine chfev ( x1, x2, f1, f2, d1, d2, ne, xe, fe, next, ierr )
!
!*******************************************************************************
!
!! CHFEV evaluates a cubic polynomial given in Hermite form.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  purpose  evaluate a cubic polynomial given in hermite form at an
!            array of points.  while designed for use by PCHFE, it may
!            be useful directly as an evaluator for a piecewise cubic
!            hermite function in applications, such as graphing, where
!            the interval is known in advance.
!
!  Description:
!
!     CHFEV evaluates the cubic polynomial determined by function values
!     f1,f2 and derivatives d1,d2 on interval (x1,x2) at the points
!     xe(j), j=1(1)ne.
!
!  Parameters:
!
!     x1,x2 -- (input) endpoints of interval of definition of cubic.
!           (error return if  x1==x2 .)
!
!     f1,f2 -- (input) values of function at x1 and x2, respectively.
!
!     d1,d2 -- (input) values of derivative at x1 and x2, respectively.
!
!     ne -- (input) number of evaluation points.  (error return if
!           ne<1 .)
!
!     xe -- (input) real array of points at which the function is to be
!           evaluated.  if any of the xe are outside the interval
!           [x1,x2], a warning error is returned in next.
!
!     fe -- (output) real array of values of the cubic function defined
!           by  x1,x2, f1,f2, d1,d2  at the points  xe.
!
!     next -- (output) integer array indicating number of extrapolation
!           points:
!            next(1) = number of evaluation points to left of interval.
!            next(2) = number of evaluation points to right of interval.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           "recoverable" errors:
!              ierr = -1  if ne<1 .
!              ierr = -2  if x1==x2 .
!                (the fe-array has not been changed in either case.)
!
  integer ne
!
  real c2
  real c3
  real d1
  real d2
  real del1
  real del2
  real delta
  real f1
  real f2
  real fe(ne)
  real h
  integer i
  integer ierr
  integer next(2)
  real x
  real x1
  real x2
  real xe(ne)
  real xma
  real xmi
!
  if (ne < 1)  go to 5001
  h = x2 - x1
  if (h == 0.0)  go to 5002
!
!  initialize.
!
  ierr = 0
  next(1) = 0
  next(2) = 0
  xmi = min ( 0.0, h )
  xma = max ( 0.0, h )
!
!  compute cubic coefficients (expanded about x1).
!
  delta = (f2 - f1)/h
  del1 = (d1 - delta)/h
  del2 = (d2 - delta)/h
!                                           (delta is no longer needed.)
  c2 = -(del1+del1 + del2)
  c3 = (del1 + del2)/h
!                               (h, del1 and del2 are no longer needed.)
!
!  evaluation loop.
!
  do i = 1, ne
     x = xe(i) - x1
     fe(i) = f1 + x*(d1 + x*(c2 + x*c3))
!          count extrapolation points.
     if ( x<xmi )  next(1) = next(1) + 1
     if ( x>xma )  next(2) = next(2) + 1
!        (note redundancy--if either condition is true, other is false.)
  end do

  return
!
!  error returns.
!
 5001 continue
!     ne<1 return.
  ierr = -1
  call xerror ('chfev -- number of evaluation points less than one', 50, &
    ierr, 1)
  return

 5002 continue
!     x1==x2 return.
  ierr = -2
  call xerror ('chfev -- interval endpoints equal', 33, ierr, 1)
  return
end
function chfiv ( x1, x2, f1, f2, d1, d2, a, b, ierr )
!
!*******************************************************************************
!
!! CHFIV evaluates the integral of a cubic polynomial in Hermite form.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!          chfiv:  cubic hermite function integral evaluator.
!
!     called by  PCHIA  to evaluate the integral of a single cubic (in
!     hermite form) over an arbitrary interval (a,b).
!
!  Parameters:
!
!     value -- (output) value of the requested integral.
!
!     x1,x2 -- (input) endpoints if interval of definition of cubic.
!           (must be distinct.  error return if not.)
!
!     f1,f2 -- (input) function values at the ends of the interval.
!
!     d1,d2 -- (input) derivative values at the ends of the interval.
!
!     a,b -- (input) endpoints of interval of integration.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0 (no errors).
!           "recoverable errors":
!              ierr = -1  if x1==x2 .
!                (value has not been set in this case.)
!
  real a
  real b
  real chfiv
  real d1
  real d2
  real dterm
  real f1
  real f2
  real fterm
  real h
  integer ierr
  real phia1
  real phia2
  real phib1
  real phib2
  real psia1
  real psia2
  real psib1
  real psib2
  real ta1
  real ta2
  real tb1
  real tb2
  real ua1
  real ua2
  real ub1
  real ub2
  real x1
  real x2
!
!  validity check input.
!
  if (x1 == x2)  go to 5001
  ierr = 0
!
!  compute integral.
!
  h = x2 - x1
  ta1 = (a - x1) / h
  ta2 = (x2 - a) / h
  tb1 = (b - x1) / h
  tb2 = (x2 - b) / h

  ua1 = ta1**3
  phia1 = ua1 * ( 2.0E+00 - ta1 )
  psia1 = ua1 * ( 3.0E+00 * ta1 - 4.0E+00 )
  ua2 = ta2**3
  phia2 =  ua2 * ( 2.0E+00 - ta2)
  psia2 = -ua2 * ( 3.0E+00 * ta2 - 4.0E+00 )

  ub1 = tb1**3
  phib1 = ub1 * ( 2.0E+00 - tb1 )
  psib1 = ub1 * ( 3.0E+00 * tb1 - 4.0E+00 )
  ub2 = tb2**3
  phib2 =  ub2 * ( 2.0E+00 - tb2 )
  psib2 = -ub2 * ( 3.0E+00 * tb2 - 4.0E+00 )

  fterm =   f1*(phia2 - phib2) + f2*(phib1 - phia1)
  dterm = ( d1*(psia2 - psib2) + d2*(psib1 - psia1) )*(h/6.0)

  chfiv = 0.5 * h * ( fterm + dterm )

  return
!
!  error return.
!
 5001 continue
  ierr = -1
  call xerror ('chfiv -- x1 equal to x2', 23, ierr, 1)
  return
end
function chfmc ( d1, d2, delta )
!
!*******************************************************************************
!
!! CHFMC determines the monotonicity properties of a cubic polynomial.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!          chfmc:  cubic hermite function monotonicity checker.
!
!    called by  pchmc  to determine the monotonicity properties of the
!    cubic with boundary derivative values d1,d2 and chord slope delta.
!
!  Parameters:
!
!     d1,d2 -- (input) derivative values at the ends of an interval.
!
!     delta -- (input) data slope over that interval.
!
!     ismon -- (output) integer function value, indicating the monoto-
!           nicity of the cubic segment:
!             ismon = -1  if function is strictly decreasing;
!             ismon =  0  if function is constant;
!             ismon =  1  if function is strictly increasing;
!             ismon =  2  if function is non-monotonic;
!             ismon =  3  if unable to determine.
!
  real a
  real b
  integer chfmc
  real d1
  real d2
  real delta
  real eps
  integer ismon
  integer itrue
  real phi
!
  eps = 10.0E+00 * epsilon ( eps )
!
!  make the check.
!
  if ( delta == 0.0E+00 )  then

    if ( d1 == 0.0E+00 .and. d2 == 0.0E+00 )  then
      ismon = 0
    else
      ismon = 2
    end if

  else

     itrue = sign ( 1.0, delta)
     a = d1/delta
     b = d2/delta
     if ((a<0.0) .or. (b<0.0))  then
        ismon = 2
     else if ((a<=3.0-eps) .and. (b<=3.0-eps))  then
!           inside square (0,3)x(0,3)  implies   ok.
        ismon = itrue
     else if ( a > 4.0+eps .and. b > 4.0+eps )  then
!           outside square (0,4)x(0,4)  implies   nonmonotonic.
        ismon = 2
     else
!           must check against boundary of ellipse.
        a = a - 2.0E+00
        b = b - 2.0E+00
        phi = ((a*a + b*b) + a*b) - 3.0E+00
        if (phi < -eps)  then
           ismon = itrue
        else if (phi > eps)  then
           ismon = 2
        else
!
!  Too close to boundary to tell,
!                  in the presence of round-off errors.
!
           ismon = 3
        end if
     end if
  end if

  chfmc = ismon

  return
end
subroutine chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )
!
!*******************************************************************************
!
!! CHKDER checks the gradients of M nonlinear functions in N variables.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!   this subroutine is a companion routine to snls1,snls1e,snsq,and
!   snsqe which may be used to check the calculation of the jacobian.
!
!   subroutine chkder
!
!     this subroutine checks the gradients of m nonlinear functions
!     in n variables, evaluated at a point x, for consistency with
!     the functions themselves. the user must call ckder twice,
!     first with mode = 1 and then with mode = 2.
!
!     mode = 1. on input, x must contain the point of evaluation.
!               on output, xp is set to a neighboring point.
!
!     mode = 2. on input, fvec must contain the functions and the
!                         rows of fjac must contain the gradients
!                         of the respective functions each evaluated
!                         at x, and fvecp must contain the functions
!                         evaluated at xp.
!               on output, err contains measures of correctness of
!                          the respective gradients.
!
!     the subroutine does not perform reliably if cancellation or
!     rounding errors cause a severe loss of significance in the
!     evaluation of a function. therefore, none of the components
!     of x should be unusually small (in particular, zero) or any
!     other value which may cause loss of significance.
!
!     the subroutine statement is
!
!     subroutine chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables.
!
!       x is an input array of length n.
!
!       fvec is an array of length m. on input when mode = 2,
!         fvec must contain the functions evaluated at x.
!
!       fjac is an m by n array. on input when mode = 2,
!         the rows of fjac must contain the gradients of
!         the respective functions evaluated at x.
!
!       ldfjac is a positive integer input parameter not less than m
!         which specifies the leading dimension of the array fjac.
!
!       xp is an array of length n. on output when mode = 1,
!         xp is set to a neighboring point of x.
!
!       fvecp is an array of length m. on input when mode = 2,
!         fvecp must contain the functions evaluated at xp.
!
!       mode is an integer input variable set to 1 on the first call
!         and 2 on the second. other values of mode are equivalent
!         to mode = 1.
!
!       err is an array of length m. on output when mode = 2,
!         err contains measures of correctness of the respective
!         gradients. if there is no severe loss of significance,
!         then if err(i) is 1.0E+00 the i-th gradient is correct,
!         while if err(i) is 0.0E+00 the i-th gradient is incorrect.
!         for values of err between 0.0E+00 and 1.0, the categorization
!         is less certain. in general, a value of err(i) greater
!         than 0.5 indicates that the i-th gradient is probably
!         correct, while a value of err(i) less than 0.5 indicates
!         that the i-th gradient is probably incorrect.
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!  references  powell, m. j. d., *a hybrid method for nonlinear
!                 equations*, 'numerical methods for nonlinear
!                 algebraic equations', p. rabinowitz, editor,
!                 gordon and breach, 1970.
!
  integer ldfjac
  integer m
  integer n
!
  real eps
  real epsf
  real epslog
  real epsmch
  real err(m)
  real, parameter :: factor = 100.0E+00
  real fjac(ldfjac,n)
  real fvec(m)
  real fvecp(m)
  integer i
  integer j
  integer mode
  real r1mach
  real temp
  real x(n)
  real xp(n)
!
  epsmch = epsilon ( epsmch )
  eps = sqrt ( epsmch )

  if ( mode == 2 ) go to 20
!
!        mode = 1.
!
    do j = 1, n
      temp = eps * abs(x(j))
      if (temp == 0.0) temp = eps
      xp(j) = x(j) + temp
    end do

     go to 70
   20 continue
!
!        mode = 2.
!
     epsf = factor*epsmch
     epslog = log10(eps)
     err(1:m) = 0.0E+00

     do j = 1, n
        temp = abs(x(j))
        if (temp == 0.0) temp = 1.0E+00
        do i = 1, m
           err(i) = err(i) + temp*fjac(i,j)
        end do
     end do

     do i = 1, m
        temp = 1.0E+00
        if ( fvec(i) /= 0.0E+00 .and. fvecp(i) /= 0.0E+00 .and. &
          abs(fvecp(i)-fvec(i)) >= epsf*abs(fvec(i)) ) then
          temp = eps*abs((fvecp(i)-fvec(i))/eps-err(i)) &
                      /(abs(fvec(i)) + abs(fvecp(i)))
        end if

        err(i) = 1.0E+00

        if ( temp > epsmch .and. temp < eps ) then
          err(i) = (log10(temp) - epslog)/epslog
        end if

        if (temp >= eps) err(i) = 0.0E+00
      end do

   70 continue

  return
end
subroutine chlhsn ( nr, n, a, epsm, sx, udiag )
!
!*******************************************************************************
!
!! CHLHSN finds the l(l-transpose) [written ll+] decomposition of the perturbed
! model hessian matrix a+mu*i(where mu\0 and i is the identity matrix)
! which is safely positive definite.
!
!  if a is safely positive definite upon entry, then mu=0.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(n,n)      <--> on entry; "a" is model hessian (only lower
!                  triangular part and diagonal stored)
!                  on exit:  a contains l of ll+ decomposition of
!                  perturbed model hessian in lower triangular
!                  part and diagonal and contains hessian in upper
!                  triangular part and udiag
! epsm         --> machine epsilon
! sx(n)        --> diagonal scaling matrix for x
! udiag(n)    <--  on exit: contains diagonal of hessian
!
! internal variables
!
! tol              tolerance
! diagmn           minimum element on diagonal of a
! diagmx           maximum element on diagonal of a
! offmax           maximum off-diagonal element of a
! offrow           sum of off-diagonal elements in a row of a
! evmin            minimum eigenvalue of a
! evmax            maximum eigenvalue of a
!
!  Description:
!
! 1. if "a" has any negative diagonal elements, then choose mu>0
! such that the diagonal of a:=a+mu*i is all positive
! with the ratio of its smallest to largest element on the
! order of sqrt(epsm).
!
! 2. "a" undergoes a perturbed cholesky decomposition which
! results in an ll+ decomposition of a+d, where d is a
! non-negative diagonal matrix which is implicitly added to
! "a" during the decomposition if "a" is not positive definite.
! "a" is retained and not changed during this process by
! copying l into the upper triangular part of "a" and the
! diagonal into udiag.  then the cholesky decomposition routine
! is called.  on return, addmax contains maximum element of d.
!
! 3. if addmax=0, "a" was positive definite going into step 2
! and return is made to calling program.  otherwise,
! the minimum number sdd which must be added to the
! diagonal of a to make it safely strictly diagonally dominant
! is calculated.  since a+addmax*i and a+sdd*i are safely
! positive definite, choose mu=min(addmax,sdd) and decompose
! a+mu*i to obtain l.
!
  integer n
  integer nr
!
  real a(nr,*)
  real addmax
  real amu
  real diagmx
  real diagmn
  real epsm
  real evmax
  real evmin
  integer i
  integer j
  real offmax
  real offro
  real offrow
  real posmax
  real sdd
  real sx(n)
  real tol
  real udiag(n)
!
! scale hessian
! pre- and post- multiply "a" by inv(sx)
!
  do j=1,n
    do i=j,n
      a(i,j)=a(i,j)/(sx(i)*sx(j))
    end do
  end do
!
! step1
!
! note:  if a different tolerance is desired throughout this
! algorithm, change tolerance here:
!
  tol=sqrt(epsm)
!
  diagmx=a(1,1)
  diagmn=a(1,1)

  do i=2,n
    if ( a(i,i)<diagmn) diagmn=a(i,i)
    if ( a(i,i)>diagmx) diagmx=a(i,i)
  end do

  posmax = max(diagmx,0.0)

  if ( diagmn > posmax*tol) go to 100

    amu=tol*(posmax-diagmn)-diagmn

    if ( amu/=0.) go to 60
!
!  Find largest off-diagonal element of A.
!
      offmax=0.0E+00

      do i=2,n
        do j=1,i-1
          if ( abs(a(i,j))>offmax) offmax=abs(a(i,j))
        end do
      end do

      amu=offmax

      if ( amu == 0.0E+00 ) then
        amu=1.0E+00
      else
        amu=amu*(1.0+tol)
      end if
!
! a=a + mu*i
!
   60 continue

    do i=1,n
      a(i,i)=a(i,i)+amu
    end do

    diagmx=diagmx+amu
!
! step2
!
! copy lower triangular part of "a" to upper triangular part
! and diagonal of "a" to udiag
!
  100 continue

  do j=1,n
    udiag(j)=a(j,j)
    do i=j+1,n
      a(j,i)=a(i,j)
    end do
  end do

  call choldc ( nr, n, a, diagmx, tol, addmax )
!
! step3
!
! if addmax=0, "a" was positive definite going into step 2,
! the ll+ decomposition has been done, and we return.
! otherwise, addmax>0.  perturb "a" so that it is safely
! diagonally dominant and find ll+ decomposition
!
  if ( addmax<=0.) go to 170
!     if ( addmax>0.)
!     then
!
! restore original "a" (lower triangular part and diagonal)
!
    do j = 1, n
      a(j,j)=udiag(j)
      do i=j+1,n
        a(i,j)=a(j,i)
      end do
    end do
!
! find sdd such that a+sdd*i is safely positive definite
! note:  evmin<0 since a is not positive definite;
!
    evmin=0.0E+00
    evmax=a(1,1)

    do i=1,n

      offrow=0.0E+00

      do j=1,i-1
        offrow=offrow+abs(a(i,j))
      end do
!
!  I'm not sure about the text of this loop.
!  It seems to have been damaged earlier, and read
!
!    offrow = offro
!
      do j=i+1,n
        offrow = offrow + abs(a(j,i))
      end do

      evmin=min(evmin,a(i,i)-offrow)
      evmax=max(evmax,a(i,i)+offrow)

    end do

    sdd=tol*(evmax-evmin)-evmin
!
! perturb "a" and decompose again
!
    amu=min(sdd,addmax)

    do i=1,n
      a(i,i)=a(i,i)+amu
      udiag(i)=a(i,i)
    end do
!
! "a" now guaranteed safely positive definite
!
    call choldc ( nr, n, a, 0.0, tol, addmax )
!
! unscale hessian and cholesky decomposition matrix
!
  170 continue

  do j=1,n

    do i=j,n
      a(i,j)=sx(i)*a(i,j)
    end do

    do i=1,j-1
      a(i,j)=sx(i)*sx(j)*a(i,j)
    end do

    udiag(j)=udiag(j)*sx(j)*sx(j)

  end do

  return
end
subroutine choldc ( nr, n, a, diagmx, tol, addmax )
!
!*******************************************************************************
!
!! CHOLDC finds the perturbed l(l-transpose) [written ll+] decomposition
! of a+d, where d is a non-negative diagonal matrix added to a if
! necessary to allow the cholesky decomposition to continue.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(n,n)      <--> on entry: matrix for which to find perturbed
!                       cholesky decomposition
!                  on exit:  contains l of ll+ decomposition
!                  in lower triangular part and diagonal of "a"
! diagmx       --> maximum diagonal element of "a"
! tol          --> tolerance
! addmax      <--  maximum amount implicitly added to diagonal of "a"
!                  in forming the cholesky decomposition of a+d
! internal variables
!
! aminl    smallest element allowed on diagonal of l
! amnlsq   =aminl**2
! offmax   maximum off-diagonal element in column of a
!
!  Description:
!
! the normal cholesky decomposition is performed.  however, if at any
! point the algorithm would attempt to set l(i,i)=sqrt(temp)
! with temp < tol*diagmx, then l(i,i) is set to sqrt(tol*diagmx)
! instead.  this is equivalent to adding tol*diagmx-temp to a(i,i)
!
  integer n
  integer nr
!
  real a(nr,n)
  real addmax
  real aminl
  real amnlsq
  real diagmx
  integer i
  integer j
  integer k
  real offmax
  real sum
  real temp
  real tol
!
  addmax = 0.0E+00
  aminl = sqrt(diagmx*tol)
  amnlsq = aminl*aminl
!
! form column j of l
!
  do j = 1, n
!
! find diagonal elements of l
!
    sum = 0.0E+00
    do k = 1,j-1
      sum = sum + a(j,k)**2
    end do

    temp = a(j,j)-sum

    if ( temp<amnlsq) go to 30

      a(j,j)=sqrt(temp)
      go to 40
!
! find maximum off-diagonal element in column
   30     offmax=0.

      do i=j+1,n
        if ( abs(a(i,j))>offmax) offmax=abs(a(i,j))
      end do

      if ( offmax<=amnlsq) offmax=amnlsq
!
! add to diagonal element  to allow cholesky decomposition to continue
!
      a(j,j)=sqrt(offmax)
      addmax=max(addmax,offmax-temp)
!
! find i,j element of lower triangular matrix
!
   40   continue

    do i=j+1,n
      sum=0.0E+00
      do k=1,j-1
        sum=sum+a(i,k)*a(j,k)
      end do
      a(i,j)=(a(i,j)-sum)/a(j,j)
    end do

  end do

  return
end
subroutine cosqb ( n, x, wsave )
!
!*******************************************************************************
!
!! COSQB computes the fast cosine transform of quarter wave data.
!
!
!  Discussion:
!
!    COSQB computes a sequence from its representation in terms of a cosine
!    series with odd wave numbers.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!
!        4 * X_in(K) * cos ( ( 2 * K - 1 ) * ( I - 1 ) * PI / ( 2 * N ) )
!
!    COSQB is the unnormalized inverse of COSQF since a call of COSQB
!    followed by a call of COSQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the algorithm.  The WSAVE array must be initialized by
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N.
!
  integer n
!
  real, parameter :: tsqrt2 = 2.82842712474619E+00
  real wsave(3*n+15)
  real x(n)
  real x1
!
  if ( n < 2 ) then
    x(1) = 4.0E+00 * x(1)
  else if ( n == 2 ) then
    x1 = 4.0E+00 * ( x(1) + x(2) )
    x(2) = tsqrt2 * ( x(1) - x(2) )
    x(1) = x1
  else
    call cosqb1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqb1 ( n, x, w, xh )
!
!*******************************************************************************
!
!! COSQB1 is a lower level routine used by COSQB.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array.
!
!    Input/output, real X(N).
!    On input, the cosine series coefficients.
!    On output, the corresponding data vector.
!
!    Input, real W(N).
!
!    Input, real XH(2*N+15).
!
  integer n
!
  integer i
  integer k
  integer kc
  integer ns2
  real w(n)
  real x(n)
  real xh(2*n+15)
  real xim1
!
  ns2 = ( n + 1 ) / 2

  do i = 3, n, 2
    xim1 = x(i-1) + x(i)
    x(i) = x(i) - x(i-1)
    x(i-1) = xim1
  end do

  x(1) = x(1) + x(1)

  if ( mod ( n, 2 ) == 0 ) then
    x(n) = 2.0E+00 * x(n)
  end if

  call rfftb ( n, x, xh )

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = w(k-1) * x(kc) + w(kc-1) * x(k)
    xh(kc) = w(k-1) * x(k) - w(kc-1) * x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * ( x(ns2+1) + x(ns2+1) )
  end if

  do k = 2, ns2
    kc = n + 2 - k
    x(k) = xh(k) + xh(kc)
    x(kc) = xh(k) - xh(kc)
  end do

  x(1) = 2.0E+00 * x(1)

  return
end
subroutine cosqf ( n, x, wsave )
!
!*******************************************************************************
!
!! COSQF computes the fast cosine transform of quarter wave data.
!
!
!  Discussion:
!
!    COSQF computes the coefficients in a cosine series representation
!    with only odd wave numbers.
!
!    COSQF is the unnormalized inverse of COSQB since a call of COSQF
!    followed by a call of COSQB will multiply the input sequence X
!    by 4*N.
!
!    The array WSAVE must be initialized by calling COSQI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + sum ( 2 <= K <= N )
!
!        2 * X_in(K) * cos ( ( 2 * I - 1 ) * ( K - 1 ) * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array X.  The method is
!    more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the algorithm.  The WSAVE array must be initialized by
!    calling COSQI.  A different WSAVE array must be used for each different
!    value of N.
!
  integer n
!
  real, parameter :: sqrt2 = 1.4142135623731E+00
  real tsqx
  real wsave(3*n+15)
  real x(n)
!
  if ( n < 2 ) then

  else if ( n == 2 ) then
    tsqx = sqrt2 * x(2)
    x(2) = x(1) - tsqx
    x(1) = x(1) + tsqx
  else
    call cosqf1 ( n, x, wsave(1), wsave(n+1) )
  end if

  return
end
subroutine cosqf1 ( n, x, w, xh )
!
!*******************************************************************************
!
!! COSQF1 is a lower level routine used by COSQF.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Input/output, real X(N).
!    On input, the data to be transformed.
!    On output, the transformed data.
!
!    Input, real W(N).
!
!    Input, real XH(2*N+15).
!
  integer n
!
  integer i
  integer k
  integer kc
  integer ns2
  real w(n)
  real x(n)
  real xh(2*n+15)
  real xim1
!
  ns2 = ( n + 1 ) / 2

  do k = 2, ns2
    kc = n + 2 - k
    xh(k) = x(k) + x(kc)
    xh(kc) = x(k) - x(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  do k = 2, ns2
    kc = n+2-k
    x(k) = w(k-1) * xh(kc) + w(kc-1) * xh(k)
    x(kc) = w(k-1) * xh(k) - w(kc-1) * xh(kc)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    x(ns2+1) = w(ns2) * xh(ns2+1)
  end if

  call rfftf ( n, x, xh )

  do i = 3, n, 2
    xim1 = x(i-1) - x(i)
    x(i) = x(i-1) + x(i)
    x(i-1) = xim1
  end do

  return
end
subroutine cosqi ( n, wsave )
!
!*******************************************************************************
!
!! COSQI initializes WSAVE, used in COSQF and COSQB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The method
!    is more efficient when N is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the COSQB and COSQF algorithms.
!
  integer n
!
  real dt
  integer k
  real pimach
  real wsave(3*n+15)
!
  dt = 0.5E+00 * pimach() / real ( n )

  do k = 1, n
    wsave(k) = cos ( real ( k ) * dt )
  end do

  call rffti ( n, wsave(n+1) )

  return
end
subroutine cost ( n, x, wsave )
!
!*******************************************************************************
!
!! COST computes the discrete Fourier cosine transform of an even sequence.
!
!
!  Discussion:
!
!    COST is the unnormalized inverse of itself since a call of COST
!    followed by another call of COST will multiply the input sequence
!    X by 2*(N-1).
!
!    The array WSAVE must be initialized by calling COSTI.
!
!    The transform is defined by:
!
!      X_out(I) = X_in(1) + (-1) **(I-1) * X_in(N) + sum ( 2 <= K <= N-1 )
!
!        2 * X_in(K) * cos ( ( K - 1 ) * ( I - 1 ) * PI / ( N - 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  The
!    method is more efficient when N-1 is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(3*N+15).
!    The WSAVE array must be initialized by calling COSTI.  A different
!    array must be used for each different value of N.
!
  integer n
!
  real c1
  integer i
  integer k
  integer kc
  integer ns2
  real t1
  real t2
  real tx2
  real wsave(3*n+15)
  real x(n)
  real x1h
  real x1p3
  real xi
  real xim2
!
  ns2 = n / 2

  if ( n <= 1 ) then
    return
  end if

  if ( n == 2 ) then
    x1h = x(1) + x(2)
    x(2) = x(1) - x(2)
    x(1) = x1h
    return
  end if

  if ( n == 3 ) then
    x1p3 = x(1) + x(3)
    tx2 = x(2) + x(2)
    x(2) = x(1) - x(3)
    x(1) = x1p3 + tx2
    x(3) = x1p3 - tx2
    return
  end if

  c1 = x(1) - x(n)
  x(1) = x(1) + x(n)

  do k = 2, ns2
    kc = n + 1 - k
    t1 = x(k) + x(kc)
    t2 = x(k) - x(kc)
    c1 = c1 + wsave(kc) * t2
    t2 = wsave(k) * t2
    x(k) = t1 - t2
    x(kc) = t1 + t2
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(ns2+1) = x(ns2+1) + x(ns2+1)
  end if

  call rfftf ( n-1, x, wsave(n+1) )

  xim2 = x(2)
  x(2) = c1

  do i = 4, n, 2
    xi = x(i)
    x(i) = x(i-2) - x(i-1)
    x(i-1) = xim2
    xim2 = xi
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n) = xim2
  end if

  return
end
subroutine costi ( n, wsave )
!
!*******************************************************************************
!
!! COSTI initializes WSAVE, used in COST.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.  The
!    method is more efficient when N-1 is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, depending on N, and
!    required by the COST algorithm.
!
  integer n
!
  real dt
  integer k
  real pimach
  real wsave(3*n+15)
!
  if ( n <= 3 ) then
    return
  end if

  dt = pimach ( ) / real ( n - 1 )

  do k = 2, ( n / 2 )
    wsave(k)     = 2.0E+00 * sin ( real ( k - 1 ) * dt )
    wsave(n+1-k) = 2.0E+00 * cos ( real ( k - 1 ) * dt )
  end do

  call rffti ( n-1, wsave(n+1) )

  return
end
function csevl ( x, cs, n )
!
!*******************************************************************************
!
!! CSEVL evaluates an N term Chebyshev series.
!
!
!  Reference:
!
!    R Broucke,
!    algorithm 446, c.a.c.m.,
!    volume 16, page 254, 1973.
!
!    Fox and Parker,
!    chebyshev polynomials in numerical analysis,
!    oxford press, page 56.
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  parameters:
!
! x    value at which the series is to be evaluated.
! cs   array of n terms of a chebyshev series.  in eval-
!      uating cs, only half the first coefficient is summed.
! n    number of terms in array cs.
!
  integer n
!
  real b0
  real b1
  real b2
  real cs(n)
  real csevl
  integer i
  real x
!
  if ( n < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CSEVL - Fatal error!'
    write ( *, * ) '  Number of terms N is less than 1.'
    stop
  end if

  if ( n > 1000 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CSEVL - Fatal error!'
    write ( *, * ) '  The number of terms is more than 1000.'
    stop
  end if

  if ( x < -1.0E+00 .or. x > 1.0E+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'CSEVL - Fatal error!'
    write ( *, * ) '  The input argument X is outside the interval [-1,1].'
    stop
  end if

  b1 = 0.0E+00
  b0 = 0.0E+00

  do i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = 2.0E+00 * x * b1 - b2 + cs(i)
  end do

  csevl = 0.5 * ( b0 - b2 )

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
subroutine d1fcn
!
!*******************************************************************************
!
!! D1FCN is a dummy routine.
!
  write ( *, * ) ' '
  write ( *, * ) 'D1FCN accidentally called!'

  stop
end
function d1mach(i)
!
!*******************************************************************************
!
!! D1MACH returns double precision machine constants.
!
!
!  Assuming that the internal representation of a double precision number is
!  in base B, with T the number of base-B digits in the mantissa, and EMIN the
!  smallest possible exponent and EMAX the largest possible exponent, then
!
!    D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!    D1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!    D1MACH(3) = B**(-T), the smallest relative spacing.
!    D1MACH(4) = B**(1-T), the largest relative spacing.
!    D1MACH(5) = log10(B).
!
!  To alter this function for a particular environment, the desired set of DATA
!  statements should be activated by removing the C from column 1.  On rare
!  machines, a STATIC statement may need to be added, but probably more systems
!  prohibit than require it.
!
!  For IEEE-arithmetic machines (binary standard), one of the first two sets of
!  constants below should be appropriate.
!
!  Where possible, octal or hexadecimal constants have been used to specify the
!  constants exactly, which has in some cases required the use of EQUIVALENCED
!  integer arrays.
!
  double precision d1mach
  integer diver(4)
  double precision dmach(5)
  integer i
  integer large(4)
  integer log10(4)
  integer right(4)
  integer small(4)
!
  equivalence (dmach(1),small(1))
  equivalence (dmach(2),large(1))
  equivalence (dmach(3),right(1))
  equivalence (dmach(4),diver(1))
  equivalence (dmach(5),log10(1))
!
!  IEEE arithmetic machines, such as the ATT 3B series and
!  Motorola 68000 based machines such as the SUN 3 and ATT PC
!  7300, and the SGI Iris, in which the most significant byte is
!  stored first.
!
   data small(1),small(2) /    1048576,          0 /
   data large(1),large(2) / 2146435071,         -1 /
   data right(1),right(2) / 1017118720,          0 /
   data diver(1),diver(2) / 1018167296,          0 /
   data log10(1),log10(2) / 1070810131, 1352628735 /
!
!  IEEE arithmetic machines and 8087-based micros, such as the IBM PC,
!  ATT 6300, DEC PMAX, DEC ALPHA, NEXT, in which the most
!  significant byte is stored last.
!
!      data small(1),small(2) /          0,    1048576 /
!      data large(1),large(2) /         -1, 2146435071 /
!      data right(1),right(2) /          0, 1017118720 /
!      data diver(1),diver(2) /          0, 1018167296 /
!      data log10(1),log10(2) / 1352628735, 1070810131 /
!
!  ALLIANT FX/8 UNIX FORTRAN compiler.
!
!      data dmach(1) / 2.22507385850721D-308 /
!      data dmach(2) / 1.79769313486231D+308 /
!      data dmach(3) / 1.1101827117665D-16 /
!      data dmach(4) / 2.2203654423533D-16 /
!      data dmach(5) / 3.01029995663981E-1 /
!
!  AMDAHL machines.
!
!      data small(1),small(2) /    1048576,          0 /
!      data large(1),large(2) / 2147483647,         -1 /
!      data right(1),right(2) /  856686592,          0 /
!      data diver(1),diver(2) /  873463808,          0 /
!      data log10(1),log10(2) / 1091781651, 1352628735 /
!
!  BURROUGHS 1700 system.
!
!      data small(1) / ZC00800000 /
!      data small(2) / Z000000000 /
!
!      data large(1) / ZDFFFFFFFF /
!      data large(2) / ZFFFFFFFFF /
!
!      data right(1) / ZCC5800000 /
!      data right(2) / Z000000000 /
!
!      data diver(1) / ZCC6800000 /
!      data diver(2) / Z000000000 /
!
!      data log10(1) / ZD00E730E7 /
!      data log10(2) / ZC77800DC0 /
!
!  BURROUGHS 5700 system.
!
!      data small(1) / O1771000000000000 /
!      data small(2) / O0000000000000000 /
!
!      data large(1) / O0777777777777777 /
!      data large(2) / O0007777777777777 /
!
!      data right(1) / O1461000000000000 /
!      data right(2) / O0000000000000000 /
!
!      data diver(1) / O1451000000000000 /
!      data diver(2) / O0000000000000000 /
!
!      data log10(1) / O1157163034761674 /
!      data log10(2) / O0006677466732724 /
!
!  BURROUGHS 6700/7700 systems.
!
!      data small(1) / O1771000000000000 /
!      data small(2) / O7770000000000000 /
!
!      data large(1) / O0777777777777777 /
!      data large(2) / O7777777777777777 /
!
!      data right(1) / O1461000000000000 /
!      data right(2) / O0000000000000000 /
!
!      data diver(1) / O1451000000000000 /
!      data diver(2) / O0000000000000000 /
!
!      data log10(1) / O1157163034761674 /
!      data log10(2) / O0006677466732724 /
!
!  CDC CYBER 170/180 series using NOS
!
!      data small(1) / O"00604000000000000000" /
!      data small(2) / O"00000000000000000000" /
!
!      data large(1) / O"37767777777777777777" /
!      data large(2) / O"37167777777777777777" /
!
!      data right(1) / O"15604000000000000000" /
!      data right(2) / O"15000000000000000000" /
!
!      data diver(1) / O"15614000000000000000" /
!      data diver(2) / O"15010000000000000000" /
!
!      data log10(1) / O"17164642023241175717" /
!      data log10(2) / O"16367571421742254654" /
!
!  CDC CYBER 170/180 series using NOS/VE
!
!      data small(1) / Z"3001800000000000" /
!      data small(2) / Z"3001000000000000" /
!
!      data large(1) / Z"4FFEFFFFFFFFFFFE" /
!      data large(2) / Z"4FFE000000000000" /
!
!      data right(1) / Z"3FD2800000000000" /
!      data right(2) / Z"3FD2000000000000" /
!
!      data diver(1) / Z"3FD3800000000000" /
!      data diver(2) / Z"3FD3000000000000" /
!
!      data log10(1) / Z"3FFF9A209A84FBCF" /
!      data log10(2) / Z"3FFFF7988F8959AC" /
!
!  CDC CYBER 200 series
!
!      data small(1) / X'9000400000000000' /
!      data small(2) / X'8FD1000000000000' /
!
!      data large(1) / X'6FFF7FFFFFFFFFFF' /
!      data large(2) / X'6FD07FFFFFFFFFFF' /
!
!      data right(1) / X'FF74400000000000' /
!      data right(2) / X'FF45000000000000' /
!
!      data diver(1) / X'FF75400000000000' /
!      data diver(2) / X'FF46000000000000' /
!
!      data log10(1) / X'FFD04D104D427DE7' /
!      data log10(2) / X'FFA17DE623E2566A' /
!
!  CDC 6000/7000 series using FTN4.
!
!      data small(1) / 00564000000000000000B /
!      data small(2) / 00000000000000000000B /
!
!      data large(1) / 37757777777777777777B /
!      data large(2) / 37157777777777777774B /
!
!      data right(1) / 15624000000000000000B /
!      data right(2) / 00000000000000000000B /
!
!      data diver(1) / 15634000000000000000B /
!      data diver(2) / 00000000000000000000B /
!
!      data log10(1) / 17164642023241175717B /
!      data log10(2) / 16367571421742254654B /
!
!  CDC 6000/7000 series using FTN5.
!
!      data small(1) / O"00564000000000000000" /
!      data small(2) / O"00000000000000000000" /
!
!      data large(1) / O"37757777777777777777" /
!      data large(2) / O"37157777777777777774" /
!
!      data right(1) / O"15624000000000000000" /
!      data right(2) / O"00000000000000000000" /
!
!      data diver(1) / O"15634000000000000000" /
!      data diver(2) / O"00000000000000000000" /
!
!      data log10(1) / O"17164642023241175717" /
!      data log10(2) / O"16367571421742254654" /
!
!  CONVEX C-1
!
!      data small(1),small(2) / '00100000'X, '00000000'X /
!      data large(1),large(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
!      data right(1),right(2) / '3CC00000'X, '00000000'X /
!      data diver(1),diver(2) / '3CD00000'X, '00000000'X /
!      data log10(1),log10(2) / '3FF34413'X, '509F79FF'X /
!
!  CONVEX C-120 (native mode) with or without -R8 option
!
!      data dmach(1) / 5.562684646268007D-309 /
!      data dmach(2) / 8.988465674311577D+307 /
!      data dmach(3) / 1.110223024625157D-016 /
!      data dmach(4) / 2.220446049250313D-016 /
!      data dmach(5) / 3.010299956639812D-001 /
!
!  CONVEX C-120 (IEEE mode) with or without -R8 option
!
!      data dmach(1) / 2.225073858507202D-308 /
!      data dmach(2) / 1.797693134862315D+308 /
!      data dmach(3) / 1.110223024625157D-016 /
!      data dmach(4) / 2.220446049250313D-016 /
!      data dmach(5) / 3.010299956639812D-001 /
!
!  CRAY 1, 2, XMP and YMP.
!
!      data small(1) / 201354000000000000000B /
!      data small(2) / 000000000000000000000B /
!
!      data large(1) / 577767777777777777777B /
!      data large(2) / 000007777777777777776B /
!
!      data right(1) / 376434000000000000000B /
!      data right(2) / 000000000000000000000B /
!
!      data diver(1) / 376444000000000000000B /
!      data diver(2) / 000000000000000000000B /
!
!      data log10(1) / 377774642023241175717B /
!      data log10(2) / 000007571421742254654B /
!
!  DATA GENERAL ECLIPSE S/200
!  Note - It may be appropriate to include the line: STATIC dmach(5)
!
!      data small /20K,3*0/
!      data large /77777K,3*177777K/
!      data right /31420K,3*0/
!      data diver /32020K,3*0/
!      data log10 /40423K,42023K,50237K,74776K/
!
!  ELXSI 6400, assuming REAL*8 is the default DOUBLE PRECISION type.
!
!      data small(1), small(2) / '00100000'X,'00000000'X /
!      data large(1), large(2) / '7FEFFFFF'X,'FFFFFFFF'X /
!      data right(1), right(2) / '3CB00000'X,'00000000'X /
!      data diver(1), diver(2) / '3CC00000'X,'00000000'X /
!      data log10(1), diver(2) / '3FD34413'X,'509F79FF'X /
!
!  HARRIS 220
!
!      data small(1),small(2) / '20000000, '00000201 /
!      data large(1),large(2) / '37777777, '37777577 /
!      data right(1),right(2) / '20000000, '00000333 /
!      data diver(1),diver(2) / '20000000, '00000334 /
!      data log10(1),log10(2) / '23210115, '10237777 /
!
!  HARRIS SLASH 6 and SLASH 7
!
!      data small(1),small(2) / '20000000, '00000201 /
!      data large(1),large(2) / '37777777, '37777577 /
!      data right(1),right(2) / '20000000, '00000333 /
!      data diver(1),diver(2) / '20000000, '00000334 /
!      data log10(1),log10(2) / '23210115, '10237777 /
!
!  HONEYWELL DPS 8/70 and 600/6000 series.
!
!      data small(1),small(2) / O402400000000, O000000000000 /
!      data large(1),large(2) / O376777777777, O777777777777 /
!      data right(1),right(2) / O604400000000, O000000000000 /
!      data diver(1),diver(2) / O606400000000, O000000000000 /
!      data log10(1),log10(2) / O776464202324, O117571775714 /
!
!  HP 2100, three word double precision option with FTN4.
!
!      data small(1), small(2), small(3) / 40000B,       0,       1 /
!      data large(1), large(2), large(3) / 77777B, 177777B, 177776B /
!      data right(1), right(2), right(3) / 40000B,       0,    265B /
!      data diver(1), diver(2), diver(3) / 40000B,       0,    276B /
!      data log10(1), log10(2), log10(3) / 46420B,  46502B,  77777B /
!
!  HP 2100, four word double precision option with FTN4.
!
!      data small(1), small(2) /  40000B,       0 /
!      data small(3), small(4) /       0,       1 /
!      data large(1), large(2) /  77777B, 177777B /
!      data large(3), large(4) / 177777B, 177776B /
!      data right(1), right(2) /  40000B,       0 /
!      data right(3), right(4) /       0,    225B /
!      data diver(1), diver(2) /  40000B,       0 /
!      data diver(3), diver(4) /       0,    227B /
!      data log10(1), log10(2) /  46420B,  46502B /
!      data log10(3), log10(4) /  76747B, 176377B /
!
!  HP 9000
!
!      d1mach(1) = 2.8480954D-306
!      d1mach(2) = 1.40444776D+306
!      d1mach(3) = 2.22044605D-16
!      d1mach(4) = 4.44089210D-16
!      d1mach(5) = 3.01029996D-1
!
!      data small(1), small(2) / 00040000000B, 00000000000B /
!      data large(1), large(2) / 17737777777B, 37777777777B /
!      data right(1), right(2) / 07454000000B, 00000000000B /
!      data diver(1), diver(2) / 07460000000B, 00000000000B /
!      data log10(1), log10(2) / 07764642023B, 12047674777B /
!
!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL SYSTEMS 85/86, PERKIN ELMER 3230,
!  and PERKIN ELMER (INTERDATA) 3230.
!
!      data small(1),small(2) / Z00100000, Z00000000 /
!      data large(1),large(2) / Z7FFFFFFF, ZFFFFFFFF /
!      data right(1),right(2) / Z33100000, Z00000000 /
!      data diver(1),diver(2) / Z34100000, Z00000000 /
!      data log10(1),log10(2) / Z41134413, Z509F79FF /
!
!  IBM PC - Microsoft FORTRAN
!
!      data small(1), small(2) / #00000000, #00100000 /
!      data large(1), large(2) / #FFFFFFFF, #7FEFFFFF /
!      data right(1), right(2) / #00000000, #3CA00000 /
!      data diver(1), diver(2) / #00000000, #3CB00000 /
!      data log10(1), log10(2) / #509F79FF, #3FD34413 /
!
!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!
!      data small(1), small(2) / Z'00000000', Z'00100000' /
!      data large(1), large(2) / Z'FFFFFFFF', Z'7FEFFFFF' /
!      data right(1), right(2) / Z'00000000', Z'3CA00000' /
!      data diver(1), diver(2) / Z'00000000', Z'3CB00000' /
!      data log10(1), log10(2) / Z'509F79FF', Z'3FD34413' /
!
!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
!  constants with Y's.
!
!      data small(1),small(2) / Z'00100000', Z'00000000' /
!      data large(1),large(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
!      data right(1),right(2) / Z'33100000', Z'00000000' /
!      data diver(1),diver(2) / Z'34100000', Z'00000000' /
!      data log10(1),log10(2) / Z'41134413', Z'509F79FF' /
!
!  PDP-10 (KA processor).
!
!      data small(1),small(2) / "033400000000, "000000000000 /
!      data large(1),large(2) / "377777777777, "344777777777 /
!      data right(1),right(2) / "113400000000, "000000000000 /
!      data diver(1),diver(2) / "114400000000, "000000000000 /
!      data log10(1),log10(2) / "177464202324, "144117571776 /
!
!  PDP-10 (KI processor).
!
!      data small(1),small(2) / "000400000000, "000000000000 /
!      data large(1),large(2) / "377777777777, "377777777777 /
!      data right(1),right(2) / "103400000000, "000000000000 /
!      data diver(1),diver(2) / "104400000000, "000000000000 /
!      data log10(1),log10(2) / "177464202324, "047674776746 /
!
!  PDP-11 FORTRANS supporting 32-bit integers (integer version).
!
!      data small(1),small(2) /    8388608,           0 /
!      data large(1),large(2) / 2147483647,          -1 /
!      data right(1),right(2) /  612368384,           0 /
!      data diver(1),diver(2) /  620756992,           0 /
!      data log10(1),log10(2) / 1067065498, -2063872008 /
!
!  PDP-11 FORTRANS supporting 32-bit integers (octal version)
!
!      data small(1),small(2) / O00040000000, O00000000000 /
!      data large(1),large(2) / O17777777777, O37777777777 /
!      data right(1),right(2) / O04440000000, O00000000000 /
!      data diver(1),diver(2) / O04500000000, O00000000000 /
!      data log10(1),log10(2) / O07746420232, O20476747770 /
!
!  PDP-11 FORTRANS supporting 16-bit integers (integer version).
!
!      data small(1),small(2) /    128,      0 /
!      data small(3),small(4) /      0,      0 /
!
!      data large(1),large(2) /  32767,     -1 /
!      data large(3),large(4) /     -1,     -1 /
!
!      data right(1),right(2) /   9344,      0 /
!      data right(3),right(4) /      0,      0 /
!
!      data diver(1),diver(2) /   9472,      0 /
!      data diver(3),diver(4) /      0,      0 /
!
!      data log10(1),log10(2) /  16282,   8346 /
!      data log10(3),log10(4) / -31493, -12296 /
!
!  PDP-11 FORTRANS supporting 16-bit integers (octal version).
!
!      data small(1),small(2) / O000200, O000000 /
!      data small(3),small(4) / O000000, O000000 /
!
!      data large(1),large(2) / O077777, O177777 /
!      data large(3),large(4) / O177777, O177777 /
!
!      data right(1),right(2) / O022200, O000000 /
!      data right(3),right(4) / O000000, O000000 /
!
!      data diver(1),diver(2) / O022400, O000000 /
!      data diver(3),diver(4) / O000000, O000000 /
!
!      data log10(1),log10(2) / O037632, O020232 /
!      data log10(3),log10(4) / O102373, O147770 /
!
!  PRIME 50 series systems with 32-bit integers and 64V MODE instructions,
!  supplied by Igor Bray.
!
!      data small(1),small(2) / :10000000000, :00000100001 /
!      data large(1),large(2) / :17777777777, :37777677775 /
!      data right(1),right(2) / :10000000000, :00000000122 /
!      data diver(1),diver(2) / :10000000000, :00000000123 /
!      data log10(1),log10(2) / :11504046501, :07674600177 /
!
!  SEQUENT BALANCE 8000
!
!      data small(1),small(2) / $00000000,  $00100000 /
!      data large(1),large(2) / $FFFFFFFF,  $7FEFFFFF /
!      data right(1),right(2) / $00000000,  $3CA00000 /
!      data diver(1),diver(2) / $00000000,  $3CB00000 /
!      data log10(1),log10(2) / $509F79FF,  $3FD34413 /
!
!  SUN Microsystems UNIX F77 compiler.
!
!      data dmach(1) / 2.22507385850720D-308 /
!      data dmach(2) / 1.79769313486231D+308 /
!      data dmach(3) / 1.1101827117665D-16 /
!      data dmach(4) / 2.2203654423533D-16 /
!      data dmach(5) / 3.01029995663981D-1 /
!
!  SUN 3 (68881 or FPA)
!
!      data small(1),small(2) / X'00100000', X'00000000' /
!      data large(1),large(2) / X'7FEFFFFF', X'FFFFFFFF' /
!      data right(1),right(2) / X'3CA00000', X'00000000' /
!      data diver(1),diver(2) / X'3CB00000', X'00000000' /
!      data log10(1),log10(2) / X'3FD34413', X'509F79FF' /
!
!  UNIVAC 1100 series.
!
!      data small(1),small(2) / O000040000000, O000000000000 /
!      data large(1),large(2) / O377777777777, O777777777777 /
!      data right(1),right(2) / O170540000000, O000000000000 /
!      data diver(1),diver(2) / O170640000000, O000000000000 /
!      data log10(1),log10(2) / O177746420232, O411757177572 /
!
!  VAX/ULTRIX F77 compiler
!
!      data small(1),small(2) /        128,           0 /
!      data large(1),large(2) /     -32769,          -1 /
!      data right(1),right(2) /       9344,           0 /
!      data diver(1),diver(2) /       9472,           0 /
!      data log10(1),log10(2) /  546979738,  -805796613 /
!
!  VAX/ULTRIX F77 compiler, G floating
!
!      data small(1), small(2) /         16,           0 /
!      data large(1), large(2) /     -32769,          -1 /
!      data right(1), right(2) /      15552,           0 /
!      data diver(1), diver(2) /      15568,           0 /
!      data log10(1), log10(2) /  1142112243, 2046775455 /
!
!  VAX-11 with FORTRAN IV-PLUS compiler
!
!      data small(1),small(2) / Z00000080, Z00000000 /
!      data large(1),large(2) / ZFFFF7FFF, ZFFFFFFFF /
!      data right(1),right(2) / Z00002480, Z00000000 /
!      data diver(1),diver(2) / Z00002500, Z00000000 /
!      data log10(1),log10(2) / Z209A3F9A, ZCFF884FB /
!
!  VAX/VMS version 2.2
!
!      data small(1),small(2) /       '80'X,        '0'X /
!      data large(1),large(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
!      data right(1),right(2) /     '2480'X,        '0'X /
!      data diver(1),diver(2) /     '2500'X,        '0'X /
!      data log10(1),log10(2) / '209A3F9A'X, 'CFF884FB'X /
!
!  VAX/VMS 11/780
!
!      data small(1), small(2) / Z00000080, Z00000000 /
!      data large(1), large(2) / ZFFFF7FFF, ZFFFFFFFF /
!      data right(1), right(2) / Z00002480, Z00000000 /
!      data diver(1), diver(2) / Z00002500, Z00000000 /
!      data log10(1), log10(2) / Z209A3F9A, ZCFF884FB /
!
!  VAX/VMS 11/780 (G-FLOATING)
!
!      data small(1), small(2) / Z00000010, Z00000000 /
!      data large(1), large(2) / ZFFFF7FFF, ZFFFFFFFF /
!      data right(1), right(2) / Z00003CC0, Z00000000 /
!      data diver(1), diver(2) / Z00003CD0, Z00000000 /
!      data log10(1), log10(2) / Z44133FF3, Z79FF509F /
!
  if(i<1.or.i>5)then
    write ( *, * ) ' '
    write(*,*)'D1MACH - Fatal error!'
    write(*,*)'I is out of bounds:',i
    d1mach=0.0d0
    stop
  else
    d1mach = dmach(i)
  end if

  return
end
subroutine d2fcn
!
!*******************************************************************************
!
!! D2FCN ??
!
  write ( *, * ) ' '
  write ( *, * ) 'd2fcn accidentally called!'
  stop
end
subroutine dfault ( n, x, typsiz, fscale, method, iexp, msg, ndigit, itnlim, &
  iagflg, iahflg, ipr, dlt, gradtl, stepmx, steptl )
!
!*******************************************************************************
!
!! DFAULT sets default values for the optimization algorithm.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! n            --> dimension of problem
! x(n)         --> initial guess to solution (to compute max step size)
! typsiz(n)   <--  typical size for each component of x
! fscale      <--  estimate of scale of minimization function
! method      <--  algorithm to use to solve minimization problem
! iexp        <--  =0 if minimization function not expensive to evaluate
! msg         <--  message to inhibit certain automatic checks + output
! ndigit      <--  number of good digits in minimization function
! itnlim      <--  maximum number of allowable iterations
! iagflg      <--  =0 if analytic gradient not supplied
! iahflg      <--  =0 if analytic hessian not supplied
! ipr         <--  device to which to send output
! dlt         <--  trust region radius
! gradtl      <--  tolerance at which gradient considered close enough
!                  to zero to terminate algorithm
! stepmx      <--  value of zero to trip default maximum in optchk
! steptl      <--  tolerance at which successive iterates considered
!                  close enough to terminate algorithm
!
  integer n
!
  real dlt
  real epsm
  real fscale
  real gradtl
  integer i1mach
  integer iagflg
  integer iahflg
  integer iexp
  integer ipr
  integer itnlim
  integer method
  integer msg
  integer ndigit
  real stepmx
  real steptl
  real typsiz(n)
  real x(n)
!
  x(n)=x(n)
!
!  Typical size of x and minimization function
!
  typsiz(1:n)=1.0E+00
  fscale=1.0E+00
!
!  Tolerances
!
  dlt=-1.0E+00
  epsm = epsilon ( epsm )
  gradtl=epsm**(1.0/3.0)
  stepmx=0.0E+00
  steptl=sqrt(epsm)
!
!  Flags
!
  method=1
  iexp=1
  msg=9
  ndigit=-1
  itnlim=150
  iagflg=0
  iahflg=0
  ipr = 6

  return
end
subroutine dogdrv ( nr, n, x, f, g, a, p, xpls, fpls, fcn, sx, stepmx, &
  steptl, dlt, iretcd, mxtake, sc, wrk1, wrk2, wrk3, ipr )
!
!*******************************************************************************
!
!! DOGDRV finds a next newton iterate (xpls) by the double dogleg method.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> old iterate x[k-1]
! f            --> function value at old iterate, f(x)
! g(n)         --> gradient  at old iterate, g(x), or approximate
! a(n,n)       --> cholesky decomposition of hessian
!                  in lower triangular part and diagonal
! p(n)         --> newton step
! xpls(n)     <--  new iterate x[k]
! fpls        <--function value at new iterate, f(xpls)
! fcn          --> name of subroutine to evaluate function
! sx(n)        --> diagonal scaling matrix for x
! stepmx       --> maximum allowable step size
! steptl       --> relative step size at which successive iterates
!                  considered close enough to terminate algorithm
! dlt         <--> trust region radius
!                  [retain value between successive calls]
! iretcd      <--  return code
!                    =0 satisfactory xpls found
!                    =1 failed to find satisfactory xpls sufficiently
!                       distinct from x
! mxtake      <--  boolean flag indicating step of maximum length used
! sc(n)        --> workspace [current step]
! wrk1(n)      --> workspace (and place holding argument to tregup)
! wrk2(n)      --> workspace
! wrk3(n)      --> workspace
! ipr          --> device to which to send output
!
  integer n
  integer nr
!
  real a(nr,n)
  real cln
  real dlt
  real eta
  real f
  real fpls
  real fplsp
  logical fstdog
  real g(n)
  integer i
  integer ipr
  integer iretcd
  logical mxtake
  logical nwtake
  real p(n)
  real rnwtln
  real sc(n)
  real stepmx
  real steptl
  real sx(n)
  real tmp
  real wrk1(n)
  real wrk2(n)
  real wrk3(n)
  real x(n)
  real xpls(n)
!
  external fcn
!
  iretcd=4
  fstdog=.true.
  tmp=0.
  do i=1,n
    tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
  end do

  rnwtln=sqrt(tmp)

  100 continue
!
! find new step by double dogleg algorithm
!
  call dogstp ( nr, n, g, a, p, sx, rnwtln, dlt, nwtake, fstdog, wrk1, &
    wrk2, cln, eta, sc, ipr, stepmx )
!
! check new point and update trust region
!
  call tregup(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl,dlt, &
    iretcd,wrk3,fplsp,xpls,fpls,mxtake,ipr,2,wrk1)

  if ( iretcd<=1) return
  go to 100
  950 format(' dogdrv    initial trust region not given.', &
             '  compute cauchy step.')
  951 format(' dogdrv    alpha =',e20.13/' dogdrv    beta  =',e20.13/ &
             ' dogdrv    dlt   =',e20.13/' dogdrv    nwtake=',l1    )
  952 format(' dogdrv    current step (sc)')
  954 format('0dogdrv    rnwtln=',e20.13)
  955 format(' dogdrv       ',5(e20.13,3x))
end
subroutine dogleg ( n, r, lr, diag, qtb, delta, x, wa1, wa2 )
!
!*******************************************************************************
!
!! DOGLEG ??
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!   subroutine dogleg
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta, the
!     problem is to determine the convex combination x of the
!     gauss-newton and scaled gradient directions that minimizes
!     (a*x - b) in the least squares sense, subject to the
!     restriction that the euclidean norm of d*x be at most delta.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization of a. that is, if a = q*r, where q has
!     orthogonal columns and r is an upper triangular matrix,
!     then dogleg expects the full upper triangle of r and
!     the first n components of (q transpose)*b.
!
!  Parameters:
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an input array of length lr which must contain the upper
!         triangular matrix r stored by rows.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       x is an output array of length n which contains the desired
!         convex combination of the gauss-newton direction and the
!         scaled gradient direction.
!
!       wa1 and wa2 are work arrays of length n.
!
  integer lr
  integer n
!
  real alpha
  real bnorm
  real delta
  real diag(n)
  real epsmch
  real gnorm
  integer i
  integer j
  integer jj
  integer jp1
  integer k
  integer l
  real qnorm
  real qtb(n)
  real r(lr)
  real r1mach
  real sgnorm
  real snrm2
  real sum
  real temp
  real wa1(n)
  real wa2(n)
  real x(n)
!
  epsmch = epsilon ( epsmch )
!
!  First, calculate the gauss-newton direction.
!
  jj = (n*(n + 1))/2 + 1

  do k = 1, n

     j = n - k + 1
     jp1 = j + 1
     jj = jj - k
     l = jj + 1
     sum = 0.0E+00

     do i = jp1, n
       sum = sum + r(l) * x(i)
       l = l + 1
     end do

     temp = r(jj)

     if ( temp == 0.0E+00 ) then

       l = j
       do i = 1, j
         temp = max ( temp, abs ( r(l) ) )
         l = l + n - i
       end do
       temp = epsmch * temp
       if (temp == 0.0) temp = epsmch

     end if

     x(j) = (qtb(j) - sum)/temp

  end do
!
!  test whether the gauss-newton direction is acceptable.
!
  wa1(1:n) = 0.0E+00
  wa2(1:n) = diag(1:n) * x(1:n)

  qnorm = snrm2(n,wa2,1)

  if (qnorm <= delta) go to 140
!
!     the gauss-newton direction is not acceptable.
!     next, calculate the scaled gradient direction.
!
  l = 1
  do j = 1, n
    temp = qtb(j)
    do i = j, n
      wa1(i) = wa1(i) + r(l)*temp
      l = l + 1
    end do
    wa1(j) = wa1(j)/diag(j)
  end do
!
!     calculate the norm of the scaled gradient direction,
!     normalize, and rescale the gradient.
!
  gnorm = snrm2(n,wa1,1)
  sgnorm = 0.0E+00
  alpha = delta/qnorm
  if (gnorm == 0.0) go to 120

  do j = 1, n
     wa1(j) = (wa1(j)/gnorm)/diag(j)
  end do
!
!     calculate the point along the scaled gradient
!     at which the quadratic is minimized.
!
  l = 1
  do j = 1, n
     sum = 0.0E+00
     do i = j, n
       sum = sum + r(l)*wa1(i)
       l = l + 1
     end do
     wa2(j) = sum
  end do
  temp = snrm2(n,wa2,1)
  sgnorm = (gnorm/temp)/temp
!
!  Test whether the scaled gradient direction is acceptable.
!
  alpha = 0.0E+00
  if (sgnorm >= delta) go to 120
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
  bnorm = snrm2(n,qtb,1)
  temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
  temp = temp - (delta/qnorm)*(sgnorm/delta)**2 &
    + sqrt((temp-(delta/qnorm))**2 &
    + (1.0-(delta/qnorm)**2)*(1.0-(sgnorm/delta)**2))

  alpha = ((delta/qnorm)*(1.0E+00 - (sgnorm/delta)**2))/temp
  120 continue
!
!     form appropriate convex combination of the gauss-newton
!     direction and the scaled gradient direction.
!
  temp = (1.0E+00 - alpha)* min (sgnorm,delta)
  do j = 1, n
    x(j) = temp*wa1(j) + alpha*x(j)
  end do

  140 continue

  return
end
subroutine dogstp ( nr, n, g, a, p, sx, rnwtln, dlt, nwtake, fstdog, ssd, v, &
  cln, eta, sc, ipr, stepmx )
!
!*******************************************************************************
!
!! DOGSTP finds a new step by double dogleg algorithm.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! g(n)         --> gradient at current iterate, g(x)
! a(n,n)       --> cholesky decomposition of hessian in
!                  lower part and diagonal
! p(n)         --> newton step
! sx(n)        --> diagonal scaling matrix for x
! rnwtln       --> newton step length
! dlt         <--> trust region radius
! nwtake      <--> boolean, =.true. if newton step taken
! fstdog      <--> boolean, =.true. if on first leg of dogleg
! ssd(n)      <--> workspace [cauchy step to the minimum of the
!                  quadratic model in the scaled steepest descent
!                  direction] [retain value between successive calls]
! v(n)        <--> workspace  [retain value between successive calls]
! cln         <--> cauchy length
!                  [retain value between successive calls]
! eta              [retain value between successive calls]
! sc(n)       <--  current step
! ipr          --> device to which to send output
! stepmx       --> maximum allowable step size
!
! internal variables
!
! cln              length of cauchy step
!
  integer n
  integer nr
!
  real a(nr,*)
  real alam
  real alpha
  real beta
  real cln
  real dlt
  real dot1
  real dot2
  real eta
  logical fstdog
  real g(n)
  integer i
  integer ipr
  integer j
  logical nwtake
  real p(n)
  real rnwtln
  real sc(n)
  real ssd(n)
  real stepmx
  real sx(n)
  real tmp
  real v(n)
!
!  Can we take a newton step?
!
  if ( rnwtln <= dlt ) then

    nwtake = .true.
    sc(1:n) = p(1:n)
    dlt = rnwtln

  else
!
! newton step too long
! cauchy step is on double dogleg curve
!
    nwtake = .false.

    if ( fstdog ) then
!
!         calculate double dogleg curve (ssd)
!
      fstdog = .false.
      alpha = 0.0E+00
      do i = 1, n
        alpha = alpha + (g(i)*g(i))/(sx(i)*sx(i))
      end do

      beta = 0.0E+00
      do i = 1, n
        tmp = 0.0E+00
        do j = i, n
          tmp = tmp + (a(j,i)*g(j))/(sx(j)*sx(j))
        end do
        beta = beta+tmp*tmp
      end do

      do i = 1, n
        ssd(i) = - (alpha/beta)*g(i)/sx(i)
      end do

      cln = alpha * sqrt ( alpha ) / beta

      eta = 0.2 + ( 0.8*alpha*alpha) / ( - beta * dot_product ( g, p ) )

      do i = 1, n
        v(i) = eta*sx(i)*p(i) - ssd(i)
      end do

      if ( dlt == - 1.0E+00 ) then
        dlt = min ( cln, stepmx )
      end if

    end if
!
!  Take a partial step in the Newton direction.
!
      if ( eta * rnwtln <= dlt ) then

        sc(1:n) = ( dlt / rnwtln ) * p(1:n)
!
!  Take a step in steepest descent direction
!
      else if ( cln >= dlt ) then

        sc(1:n) = ( dlt / cln ) * ssd(1:n) / sx(1:n)
!
!  convex combination of ssd and eta*p which has scaled length dlt
!
      else

        dot1 = dot_product ( v, ssd )
        dot2 = dot_product ( v, v )
        alam = (-dot1+sqrt((dot1*dot1)-dot2*(cln*cln-dlt*dlt)))/dot2

        sc(1:n) = ( ssd(1:n) + alam * v(1:n) ) / sx(1:n)

      end if

  end if

  return
end
subroutine ea ( newflg, svalue, limexp, result, abserr, epstab, ierr )
!
!*******************************************************************************
!
!! EA performs extrapolation to accelerate the convergence of a sequence.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  purpose  given a slowly convergent sequence, this routine attempts
!            to extrapolate nonlinearly to a better estimate of the
!            sequence's limiting value, thus improving the rate of
!            convergence. routine is based on the epsilon algorithm
!            of p. wynn. an estimate of the absolute error is also
!            given.
!
  real abserr
  real delta1
  real delta2
  real delta3
  real eprn
  real epstab(*)
  real error
  real err1
  real err2
  real err3
  real e0
  real e1
  real e2
  real e3
  integer i
  integer ib
  integer ib2
  integer ie
  integer ierr
  integer in
  integer k1
  integer k2
  integer k3
  integer limexp
  integer n
  integer newelm
  integer num
  integer nres
  logical newflg
  real relpr
  real res
  real result
  real res3la(3)
  real r1mach
  real ss
  real svalue
  real tol1
  real tol2
  real tol3
!
  if ( limexp < 3 ) then
    ierr = 1
    call xerror('limexp is less than 3',21,1,1)
    go to 110
  end if

  ierr = 0
  res3la(1)=epstab(limexp+5)
  res3la(2)=epstab(limexp+6)
  res3la(3)=epstab(limexp+7)
  result=svalue

  if ( newflg) then
    n=1
    nres=0
    newflg=.false.
    epstab(n)=svalue
    abserr=abs(result)
    go to 100
  else
    n=int(epstab(limexp+3))
    nres=int(epstab(limexp+4))
    if ( n==2) then
      epstab(n)=svalue
      abserr=6.0*abs(result-epstab(1))
      go to 100
    end if
  end if

  epstab(n)=svalue
  relpr=r1mach(4)
  eprn=10.0*relpr
  epstab(n+2)=epstab(n)
  newelm=(n-1)/2
  num=n
  k1=n

  do i=1,newelm

    k2=k1-1
    k3=k1-2
    res=epstab(k1+2)
    e0=epstab(k3)
    e1=epstab(k2)
    e2=res
    delta2=e2-e1
    err2=abs(delta2)
    tol2=max(abs(e2),abs(e1))*relpr
    delta3=e1-e0
    err3=abs(delta3)
    tol3=max(abs(e1),abs(e0))*relpr

!
!  if e0, e1 and e2 are equal to within machine accuracy, convergence is assumed.
!
    if ( err2 > tol2 .or. err3 > tol3 ) go to 10
    result=res
    abserr=err2+err3
    go to 50

   10   continue

    if ( i/=1) then
      e3=epstab(k1)
      epstab(k1)=e1
      delta1=e1-e3
      err1=abs(delta1)
      tol1=max(abs(e1),abs(e3))*relpr
!
!           if two elements are very close to each other, omit
!           a part of the table by adjusting the value of n
!
      if ( err1<=tol1.or.err2<=tol2.or.err3<=tol3) go to 20
      ss=1.0/delta1+1.0/delta2-1.0/delta3
    else
      epstab(k1)=e1
      if ( err2<=tol2.or.err3<=tol3) go to 20
      ss=1.0/delta2-1.0/delta3
    end if
!
!           test to detect irregular behaviour in the table, and
!           eventually omit a part of the table adjusting the value
!           of n
!
    if ( abs(ss*e1)>0.1e-03) go to 30
   20   n=i+i-1

    if ( nres==0) then
      abserr=err2+err3
      result=res
    else if ( nres==1) then
      result=res3la(1)
    else if ( nres==2) then
      result=res3la(2)
    else
      result=res3la(3)
    end if

    go to 50
!
!  compute a new element and eventually adjust the value of result
!
   30   res=e1+1.0/ss
    epstab(k1)=res
    k1=k1-2
    if ( nres==0) then
      abserr=err2+abs(res-e2)+err3
      result=res
      go to 40
    else if ( nres==1) then
      error=6.0*(abs(res-res3la(1)))
    else if ( nres==2) then
      error=2.0*(abs(res-res3la(2))+abs(res-res3la(1)))
    else
      error=abs(res-res3la(3))+abs(res-res3la(2))+abs(res-res3la(1))
    end if

    if ( error <= 10.0E+00 * abserr ) then
      abserr = error
      result = res
    end if

   40 continue

  end do
!
!  compute error estimate
!
    if ( nres==1) then
      abserr=6.0*(abs(result-res3la(1)))
    else if ( nres==2) then
      abserr=2.0*abs(result-res3la(2))+abs(result-res3la(1))
    else if ( nres>2) then
      abserr=abs(result-res3la(3))+abs(result-res3la(2)) &
        +abs(result-res3la(1))
    end if
!
!  shift the table
!
   50 if ( n==limexp) n=2*(limexp/2)-1
  ib=1
  if ( (num/2)*2==num) ib=2
  ie=newelm+1

  do i = 1, ie
    ib2=ib+2
    epstab(ib)=epstab(ib2)
    ib=ib2
  end do

  if ( num==n) go to 80
  in=num-n+1

  do i=1,n
    epstab(i)=epstab(in)
    in=in+1
  end do
!
!  update res3la
!
   80 continue

  if ( nres==0) then
    res3la(1)=result
  else if ( nres==1) then
    res3la(2)=result
  else if ( nres==2) then
    res3la(3)=result
  else
    res3la(1)=res3la(2)
    res3la(2)=res3la(3)
    res3la(3)=result
  end if

90 continue

  abserr=max(abserr,eprn*abs(result))
  nres=nres+1

100 continue

  n=n+1
  epstab(limexp+3)=real(n)
  epstab(limexp+4)=real(nres)
  epstab(limexp+5)=res3la(1)
  epstab(limexp+6)=res3la(2)
  epstab(limexp+7)=res3la(3)

110 continue

  return
end
function erf ( x )
!
!*******************************************************************************
!
!! ERF computes the error function.
!
!
!  Definition:
!
!    ERF(X) = ( 2 / SQRT ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( -T**2 ) dT
!
!  Modified:
!
!    10 November 1999
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, real X, the argument of the error function.
!
!    Output, real ERF, the value of the error function at X.
!
  real, parameter :: sqrtpi = 1.7724538509055160
!
  real csevl
  real erf
  real erfc
  real erfcs(13)
  integer inits
  integer, save :: nterf = 0
  real r1mach
  real sqeps
  real x
  real, save :: xbig = 0.0E+00
  real y
!
  data erfcs( 1) /   -0.049046121234691808 /
  data erfcs( 2) /   -0.14226120510371364 /
  data erfcs( 3) /    0.010035582187599796 /
  data erfcs( 4) /   -0.000576876469976748 /
  data erfcs( 5) /    0.000027419931252196 /
  data erfcs( 6) /   -0.000001104317550734 /
  data erfcs( 7) /    0.000000038488755420 /
  data erfcs( 8) /   -0.000000001180858253 /
  data erfcs( 9) /    0.000000000032334215 /
  data erfcs(10) /   -0.000000000000799101 /
  data erfcs(11) /    0.000000000000017990 /
  data erfcs(12) /   -0.000000000000000371 /
  data erfcs(13) /    0.000000000000000007 /
  data sqeps / 0.0E+00 /
!
  if ( nterf == 0 ) then
    nterf = inits ( erfcs, 13, 0.1*r1mach(3) )
    xbig = sqrt ( - log ( sqrtpi * r1mach(3) ) )
    sqeps = sqrt ( 2.0E+00 * epsilon ( sqeps ) )
  end if

  y = abs ( x )

  if ( y <= sqeps ) then
    erf = 2.0E+00 * x / sqrtpi
  else if ( y <= 1.0E+00 ) then
    erf = x * ( 1.0E+00 + csevl ( 2.0E+00 * x**2 - 1.0, erfcs, nterf ) )
  else if ( y <= xbig ) then
    erf = sign ( 1.0E+00 - erfc ( y ), x )
  else
    erf = sign ( 1.0, x )
  end if

  return
end
function erfc ( x )
!
!*******************************************************************************
!
!! ERFC computes the complementary error function.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! erfc(x) calculates the single precision complementary error
! function for single precision argument x.
!
! series for erf        on the interval  0.          to  1.00000d+00
!                                        with weighted error   7.10e-18
!                                         log weighted error  17.15
!                               significant figures required  16.31
!                                    decimal places required  17.71
!
! series for erfc       on the interval  0.          to  2.50000d-01
!                                        with weighted error   4.81e-17
!                                         log weighted error  16.32
!                        approx significant figures required  15.0E+00
!
!
! series for erc2       on the interval  2.50000d-01 to  1.00000d+00
!                                        with weighted error   5.22e-17
!                                         log weighted error  16.28
!                        approx significant figures required  15.0E+00
!                                    decimal places required  16.96
!
  real csevl
  real erfc
  real erfcs(13)
  real erfccs(24)
  real erc2cs(23)
  real eta
  integer inits
  integer nterc2
  integer nterf
  integer nterfc
  real r1mach
  real sqeps
  real sqrtpi
  real x
  real xmax
  real xsml
  real y
!
  data erfcs( 1) /   -.049046121234691808e0 /
  data erfcs( 2) /   -.14226120510371364e0 /
  data erfcs( 3) /    .010035582187599796e0 /
  data erfcs( 4) /   -.000576876469976748e0 /
  data erfcs( 5) /    .000027419931252196e0 /
  data erfcs( 6) /   -.000001104317550734e0 /
  data erfcs( 7) /    .000000038488755420e0 /
  data erfcs( 8) /   -.000000001180858253e0 /
  data erfcs( 9) /    .000000000032334215e0 /
  data erfcs(10) /   -.000000000000799101e0 /
  data erfcs(11) /    .000000000000017990e0 /
  data erfcs(12) /   -.000000000000000371e0 /
  data erfcs(13) /    .000000000000000007e0 /
  data erc2cs( 1) /   -.069601346602309501e0 /
  data erc2cs( 2) /   -.041101339362620893e0 /
  data erc2cs( 3) /    .003914495866689626e0 /
  data erc2cs( 4) /   -.000490639565054897e0 /
  data erc2cs( 5) /    .000071574790013770e0 /
  data erc2cs( 6) /   -.000011530716341312e0 /
  data erc2cs( 7) /    .000001994670590201e0 /
  data erc2cs( 8) /   -.000000364266647159e0 /
  data erc2cs( 9) /    .000000069443726100e0 /
  data erc2cs(10) /   -.000000013712209021e0 /
  data erc2cs(11) /    .000000002788389661e0 /
  data erc2cs(12) /   -.000000000581416472e0 /
  data erc2cs(13) /    .000000000123892049e0 /
  data erc2cs(14) /   -.000000000026906391e0 /
  data erc2cs(15) /    .000000000005942614e0 /
  data erc2cs(16) /   -.000000000001332386e0 /
  data erc2cs(17) /    .000000000000302804e0 /
  data erc2cs(18) /   -.000000000000069666e0 /
  data erc2cs(19) /    .000000000000016208e0 /
  data erc2cs(20) /   -.000000000000003809e0 /
  data erc2cs(21) /    .000000000000000904e0 /
  data erc2cs(22) /   -.000000000000000216e0 /
  data erc2cs(23) /    .000000000000000052e0 /
  data erfccs( 1) /   0.0715179310202925e0 /
  data erfccs( 2) /   -.026532434337606719e0 /
  data erfccs( 3) /    .001711153977920853e0 /
  data erfccs( 4) /   -.000163751663458512e0 /
  data erfccs( 5) /    .000019871293500549e0 /
  data erfccs( 6) /   -.000002843712412769e0 /
  data erfccs( 7) /    .000000460616130901e0 /
  data erfccs( 8) /   -.000000082277530261e0 /
  data erfccs( 9) /    .000000015921418724e0 /
  data erfccs(10) /   -.000000003295071356e0 /
  data erfccs(11) /    .000000000722343973e0 /
  data erfccs(12) /   -.000000000166485584e0 /
  data erfccs(13) /    .000000000040103931e0 /
  data erfccs(14) /   -.000000000010048164e0 /
  data erfccs(15) /    .000000000002608272e0 /
  data erfccs(16) /   -.000000000000699105e0 /
  data erfccs(17) /    .000000000000192946e0 /
  data erfccs(18) /   -.000000000000054704e0 /
  data erfccs(19) /    .000000000000015901e0 /
  data erfccs(20) /   -.000000000000004729e0 /
  data erfccs(21) /    .000000000000001432e0 /
  data erfccs(22) /   -.000000000000000439e0 /
  data erfccs(23) /    .000000000000000138e0 /
  data erfccs(24) /   -.000000000000000048e0 /
  data sqrtpi /1.7724538509055160e0/
  data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0./
!
  if ( nterf == 0 ) then

    eta = 0.1 * r1mach(3)
    nterf = inits ( erfcs, 13, eta )
    nterfc = inits ( erfccs, 24, eta )
    nterc2 = inits ( erc2cs, 23, eta )

    xsml = -sqrt ( - log ( sqrtpi * r1mach(3) ) )
    xmax = sqrt ( - log ( sqrtpi * r1mach(1) ) )
    xmax = xmax - 0.5 * log ( xmax ) / xmax - 0.01
    sqeps = sqrt ( 2.0E+00 * r1mach(3) )

  end if

  if ( x <= xsml ) then
    erfc = 2.0E+00
    return
  end if

  if ( x > xmax ) go to 40
  y = abs(x)
  if (y>1.0) go to 30
!
!  erfc(x) = 1.0E+00 - erf(x) for -1. <= x <= 1.
!
  if ( y < sqeps ) then
    erfc = 1.0E+00 - 2.0*x/sqrtpi
  else if (y>=sqeps) then
    erfc = 1.0E+00 - x*(1.0E+00 + csevl (2.*x*x-1., erfcs, nterf) )
  end if

  return
!
! erfc(x) = 1.0E+00 - erf(x) for 1. < abs(x) <= xmax
!
30 continue

  y = y*y

  if (y<=4.) then
    erfc = exp(-y)/abs(x) * (0.5 + csevl ((8./y-5.)/3.,erc2cs, nterc2) )
  else
    erfc = exp(-y)/abs(x) * (0.5 + csevl (8./y-1.,erfccs, nterfc) )
  end if

  if ( x < 0.0E+00 ) erfc = 2.0E+00 - erfc

  return

40  continue

  call xerror ( 'erfc    x so big erfc underflows', 32, 1, 1)
  erfc = 0.0E+00

  return
end
subroutine ezfftb ( n, r, azero, a, b, wsave )
!
!*******************************************************************************
!
!! EZFFTB computes a real periodic sequence from its Fourier coefficients.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    EZFFTB is a simplified but slower version of RFFTB.
!
!    The transform is defined by:
!
!      R(I) = AZERO + sum ( 1 <= K <= N/2 )
!
!          A(K) * cos ( K * ( I - 1 ) * 2 * PI / N )
!        + B(K) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the output array.  The
!    method is more efficient when N is the product of small primes.
!
!    Output, real R(N), the reconstructed data sequence.
!
!    Input, real AZERO, the constant Fourier coefficient.
!
!    Input, real A(N/2), B(N/2), the Fourier coefficients.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling EZFFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
  integer n
!
  real a(n/2)
  real azero
  real b(n/2)
  integer i
  integer ns2
  real r(n)
  real wsave(3*n+15)
!
  if ( n < 2 ) then

    r(1) = azero

  else if ( n == 2 ) then

    r(1) = azero + a(1)
    r(2) = azero - a(1)

  else

    ns2 = ( n - 1 ) / 2

    do i = 1, ns2
      r(2*i) = 0.5E+00 * a(i)
      r(2*i+1) = -0.5E+00 * b(i)
    end do

    r(1) = azero

    if ( mod ( n, 2 ) == 0 ) then
      r(n) = a(ns2+1)
    end if

    call rfftb ( n, r, wsave(n+1) )

  end if

  return
end
subroutine ezfftf ( n, r, azero, a, b, wsave )
!
!*******************************************************************************
!
!! EZFFTF computes the Fourier coefficients of a real periodic sequence.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    EZFFTF is a simplified but slower version of RFFTF.
!
!    The transform is defined by:
!
!      AZERO = sum ( 1 <= I <= N ) R(I) / N,
!
!    and, for K = 1 to (N-1)/2,
!
!      A(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * cos ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      A(N/2) = sum ( 1 <= I <= N ) (-1) **(I-1) * R(I) / N
!
!    For K = 1 to (N-1)/2,
!
!      B(K) = sum ( 1 <= I <= N )
!        ( 2 / N ) * R(I) * sin ( K * ( I - 1 ) * 2 * PI / N )
!
!    and, if N is even, then
!
!      B(N/2) = 0.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input, real R(N), the sequence to be transformed.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling EZFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
!    Output, real AZERO, the constant Fourier coefficient.
!
!    Output, real A(N/2), B(N/2), the Fourier coefficients.
!
  integer n
!
  real a(n/2)
  real azero
  real b(n/2)
  real cf
  integer i
  integer ns2
  real r(n)
  real wsave(3*n+15)
!
  if ( n < 2 ) then

    azero = r(1)

  else if ( n == 2 ) then

    azero = 0.5E+00 * ( r(1) + r(2) )
    a(1) = 0.5E+00 * ( r(1) - r(2) )

  else

    wsave(1:n) = r(1:n)

    call rfftf ( n, wsave(1), wsave(n+1) )

    cf = 2.0E+00 / real ( n )
    azero = 0.5E+00 * cf * wsave(1)
    ns2 = ( n + 1 ) / 2

    do i = 1, ns2-1
      a(i) = cf * wsave(2*i)
      b(i) = -cf * wsave(2*i+1)
    end do

    if ( mod ( n, 2 ) /= 1 ) then
      a(ns2) = 0.5E+00 * cf * wsave(n)
      b(ns2) = 0.0E+00
    end if

  end if

  return
end
subroutine ezffti ( n, wsave )
!
!*******************************************************************************
!
!! EZFFTI initializes WSAVE, used in EZFFTF and EZFFTB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Output, real WSAVE(3*N+15), contains data, dependent on the value
!    of N, which is necessary for the EZFFTF or EZFFTB routines.
!
  integer n
!
  real wsave(3*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call ezffti1 ( n, wsave(2*n+1), wsave(3*n+1) )

  return
end
subroutine ezffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! EZFFTI1 is a lower level routine used by EZFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Output, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  real arg1
  real argh
  real ch1
  real ch1h
  real dch1
  real dsh1
  integer i
  integer ib
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer is
  integer j
  integer k1
  integer l1
  integer l2
  integer nf
  real pimach
  real sh1
  real wa(n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0E+00 * pimach() / real ( n )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    l2 = l1 * ip
    ido = n / l2
    arg1 = real ( l1 ) * argh
    ch1 = 1.0E+00
    sh1 = 0.0E+00
    dch1 = cos ( arg1 )
    dsh1 = sin ( arg1 )

    do j = 1, ip-1

      ch1h = dch1 * ch1 - dsh1 * sh1
      sh1  = dch1 * sh1 + dsh1 * ch1
      ch1 = ch1h
      i = is + 2
      wa(i-1) = ch1
      wa(i) = sh1

      do ii = 5, ido, 2
        i = i + 2
        wa(i-1) = ch1 * wa(i-3) - sh1 * wa(i-2)
        wa(i)   = ch1 * wa(i-2) + sh1 * wa(i-3)
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
subroutine fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )
!
!*******************************************************************************
!
!! FDJAC1 computes a forward-difference approximation to a jacobian matrix.
!
!
!  Discussion:
!
!    The n by n jacobian matrix is associated with a specified
!    problem of n functions in n variables. if the jacobian has
!    a banded form, then function evaluations are saved by only
!    approximating the nonzero terms.
!
!  Modified:
!
!    25 August 2000
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         real x(n),fvec(n)
!
!         calculate the functions at x and
!         return this vector in fvec.
!
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of FDJAC1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an input array of length n.
!
!       fvec is an input array of length n which must contain the
!       functions evaluated at x.
!
!       fjac is an output n by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of FJDAC1. see description of fcn.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!       functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
  integer ldfjac
  integer n
!
  real eps
  real epsfcn
  real epsmch
  real fjac(ldfjac,n)
  real fvec(n)
  real h
  integer i
  integer iflag
  integer j
  integer k
  integer ml
  integer msum
  integer mu
  real r1mach
  real temp
  real wa1(n)
  real wa2(n)
  real x(n)
!
  external fcn
!
  epsmch = epsilon ( 1.0E+00 )
!
  eps = sqrt( max (epsfcn,epsmch))
  msum = ml + mu + 1
  if (msum < n) go to 40
!
!  computation of dense approximate jacobian.
!
     do j = 1, n

       temp = x(j)
       h = eps*abs(temp)
       if (h == 0.0) h = eps
       x(j) = temp + h
       call fcn(n,x,wa1,iflag)
       if (iflag < 0) go to 30
       x(j) = temp
       fjac(1:n,j) = (wa1(1:n) - fvec(1:n))/h

    end do

   30    continue

     go to 110

   40 continue
!
!  Computation of banded approximate jacobian.
!
     do k = 1, msum

        do j = k, n, msum
           wa2(j) = x(j)
           h = eps*abs(wa2(j))
           if (h == 0.0) h = eps
           x(j) = wa2(j) + h
        end do

        call fcn(n,x,wa1,iflag)
        if (iflag < 0) go to 100

        do j = k, n, msum

           x(j) = wa2(j)
           h = eps*abs(wa2(j))
           if (h == 0.0) h = eps

           do i = 1, n
              fjac(i,j) = 0.0E+00
              if (i >= j - mu .and. i <= j + ml) then
                fjac(i,j) = (wa1(i) - fvec(i))/h
              end if
           end do

        end do

      end do

  100    continue
  110 continue

  return
end
function fmin ( ax, bx, f, tol )
!
!*******************************************************************************
!
!! FMIN seeks a minimizer of a scalar function of a scalar variable.
!
!
!  Discussion:
!
!    FMIN seeks an approximation to the point where F attains a minimum on
!    the interval (ax,bx).
!
!     the method used is a combination of golden section search and
!     successive parabolic interpolation.  convergence is never much
!     slower than that for a fibonacci search.  if f has a continuous
!     second derivative which is positive at the minimum (which is not
!     at ax or bx), then convergence is superlinear, and usually of the
!     order of about 1.324....
!
!     the function f is never evaluated at two points closer together
!     than eps*abs(fmin) + (tol/3), where eps is approximately the
!     square root of the relative machine precision.  if f is a unimodal
!     function and the computed values of f are always unimodal when
!     separated by at least eps*abs(xstar) + (tol/3), then fmin
!     approximates the abcissa of the global minimum of f on the
!     interval ax,bx with an error less than 3*eps*abs(fmin) + tol.
!     if f is not unimodal, then fmin may approximate a local, but
!     perhaps non-global, minimum to the same accuracy.
!
!     this function subprogram is a slightly modified version of the
!     algol 60 procedure localmin given in richard brent, algorithms for
!     minimization without derivatives, prentice-hall, inc. (1973).
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization without Derivatives,
!    Prentice Hall, 1973.
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! input parameters
!
!  ax    (real)  left endpoint of initial interval
!  bx    (real) right endpoint of initial interval
!  f     real function of the form real function f(x) which evaluates
!          f(x)  for any  x in the interval  (ax,bx)
!        must be declared external in calling routine.
!  tol   (real) desired length of the interval of uncertainty of the
!        final result ( >= 0.0)
!
!
! output parameters
!
! fmin   abcissa approximating the minimizer of f.
!
! ax     lower bound for minimizer.
!
! bx     upper bound for minimizer
!
  real a
  real ax
  real b
  real bx
  real c
  real d
  real e
  real eps
  real f
  real fmin
  real fu
  real fv
  real fw
  real fx
  real p
  real q
  real r
  real tol
  real tol1
  real tol2
  real u
  real v
  real w
  real x
  real xm
!
  external f
!
  c = 0.5 * ( 3.0E+00 - sqrt ( 5.0E+00 ) )
!
!  C is the squared inverse of the golden ratio.
!
!  EPS is the square root of the relative machine precision.
!
  eps = sqrt ( epsilon ( eps ) )
!
!  initialization
!
  a = ax
  b = bx
  v = a + c*(b - a)
  w = v
  x = v
  e = 0.0E+00
  fx = f(x)
  fv = fx
  fw = fx
!
!  main loop starts here
!
   20 continue

  xm = 0.5 * ( a + b )
  tol1 = eps*abs(x) + tol/3.0E+00
  tol2 = 2.0*tol1
!
!  Check the stopping criterion.
!
  if (abs(x - xm) <= (tol2 - 0.5*(b - a))) go to 90
!
! is golden-section necessary
!
  if (abs(e) <= tol1) go to 40
!
!  fit parabola
!
  r = (x - w)*(fx - fv)
  q = (x - v)*(fx - fw)
  p = (x - v)*q - (x - w)*r
  q = 2.0*(q - r)
  if (q > 0.0) p = -p
  q = abs(q)
  r = e
  e = d
!
!  is parabola acceptable
!
   30 continue

  if (abs(p) >= abs(0.5*q*r)) go to 40
  if (p <= q*(a - x)) go to 40
  if (p >= q*(b - x)) go to 40
!
!  a parabolic interpolation step
!
  d = p/q
  u = x + d
!
!  f must not be evaluated too close to ax or bx
!
  if ((u - a) < tol2) d = sign(tol1, xm - x)
  if ((b - u) < tol2) d = sign(tol1, xm - x)
  go to 50
!
!  a golden-section step
!
   40 continue

  if (x >= xm) e = a - x
  if (x < xm) e = b - x
  d = c*e
!
!  f must not be evaluated too close to x
!
   50 if (abs(d) >= tol1) u = x + d
  if (abs(d) < tol1) u = x + sign(tol1, d)
  fu = f(u)
!
!  update  a, b, v, w, and x
!
  if (fu > fx) go to 60

  if (u >= x) a = x
  if (u < x) b = x
  v = w
  fv = fw
  w = x
  fw = fx
  x = u
  fx = fu
  go to 20

60 continue

  if (u < x) a = u
  if (u >= x) b = u
  if (fu <= fw) go to 70
  if (w == x) go to 70
  if (fu <= fv) go to 80
  if (v == x) go to 80
  if (v == w) go to 80
  go to 20

   70 v = w
  fv = fw
  w = u
  fw = fu
  go to 20
   80 v = u
  fv = fu
  go to 20
!
!end of main loop
!
   90 fmin = x
  return
end
subroutine forslv ( nr, n, a, x, b )
!
!*******************************************************************************
!
!! FORSLV solves A*x=b where A is lower triangular matrix.
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(n,n)       --> lower triangular matrix (preserved)
! x(n)        <--  solution vector
! b(n)         --> right-hand side vector
!
! note
!
! if b is no longer required by calling routine,
! then vectors b and x may share the same storage.
!
  integer n
  integer nr
!
  real a(nr,n)
  real b(n)
  integer i
  integer j
  real sum
  real x(n)
!
  x(1) = b(1) / a(1,1)

  do i = 2, n
    sum = 0.0E+00
    do j = 1, i-1
      sum = sum + a(i,j) * x(j)
    end do
    x(i) = ( b(i) - sum ) / a(i,i)
  end do

  return
end
subroutine fstocd ( n, x, fcn, sx, rnoise, g )
!
!*******************************************************************************
!
!! FSTOCD approximates the gradient of a function using central differences.
!
!
! parameters
!
! n            --> dimension of problem
! x            --> point at which gradient is to be approximated.
! fcn          --> name of subroutine to evaluate function.
! sx           --> diagonal scaling matrix for x.
! rnoise       --> relative noise in fcn [f(x)].
! g           <--  central difference approximation to gradient.
!
!
  integer n
!
  real fminus
  real fplus
  real g(n)
  integer i
  real rnoise
  real stepi
  real sx(n)
  real third
  real x(n)
  real xtempi
!
  external fcn
!
! find i th  stepsize, evaluate two neighbors in direction of i th
! unit vector, and evaluate i th  component of gradient.
!
  third = 1.0 /3.0E+00

  do i = 1, n
     stepi = rnoise**third * max(abs(x(i)), 1.0 /sx(i))
     xtempi = x(i)
     x(i) = xtempi + stepi
     call fcn (n, x, fplus)
     x(i) = xtempi - stepi
     call fcn (n, x, fminus)
     x(i) = xtempi
     g(i) = (fplus - fminus)/(2.0 *stepi)
  end do

  return
end
subroutine fstofd ( nr, m, n, xpls, fcn, fpls, a, sx, rnoise, fhat, icase )
!
!*******************************************************************************
!
!! FSTOFD finds first order forward finite difference approximation "a" to the
! first derivative of the function defined by the subprogram "fname"
! evaluated at the new iterate "xpls".
!
!
! for optimization use this routine to estimate:
! 1) the first derivative (gradient) of the optimization function "fcn
!    analytic user routine has been supplied;
! 2) the second derivative (hessian) of the optimization function
!    if no analytic user routine has been supplied for the hessian but
!    one has been supplied for the gradient ("fcn") and if the
!    optimization function is inexpensive to evaluate
!
! note
!
! _m=1 (optimization) algorithm estimates the gradient of the function
!      (fcn).   fcn(x) # f: r(n)-->r(1)
! _m=n (systems) algorithm estimates the jacobian of the function
!      fcn(x) # f: r(n)-->r(n).
! _m=n (optimization) algorithm estimates the hessian of the optimizatio
!    function, where the hessian is the first derivative of "fcn"
!
!  Parameters:
!
! nr           --> row dimension of matrix
! m            --> number of rows in a
! n            --> number of columns in a; dimension of problem
! xpls(n)      --> new iterate:  x[k]
! fcn          --> name of subroutine to evaluate function
! fpls(m)      --> _m=1 (optimization) function value at new iterate:
!                       fcn(xpls)
!                  _m=n (optimization) value of first derivative
!                       (gradient) given by user function fcn
!                  _m=n (systems)function value of associated
!                       minimization function
! a(nr,n)     <--  finite difference approximation (see note).  only
!                  lower triangular matrix and diagonal are returned
! sx(n)        --> diagonal scaling matrix for x
! rnoise       --> relative noise in fcn [f(x)]
! fhat(m)      --> workspace
! icase        --> =1 optimization (gradient)
!                  =2 systems
!                  =3 optimization (hessian)
!
! internal variables
!
! stepsz - stepsize in the j-th variable direction
!
  integer m
  integer n
  integer nr
!
  real a(nr,n)
  real fhat(m)
  real fpls(m)
  integer i
  integer icase
  integer j
  real rnoise
  real stepsz
  real sx(n)
  real xpls(n)
  real xtmpj
!
  external fcn
!
! find j-th column of a
! each column is derivative of f(fcn) with respect to xpls(j)
!
  do j=1,n
    stepsz=sqrt(rnoise)*max(abs(xpls(j)),1./sx(j))
    xtmpj=xpls(j)
    xpls(j)=xtmpj+stepsz
    call fcn(n,xpls,fhat)
    xpls(j)=xtmpj
    do i=1,m
      a(i,j)=(fhat(i)-fpls(i))/stepsz
    end do
  end do

  if ( icase/=3) then
    return
  end if
!
! if computing hessian, a must be symmetric
!
  do j=1,n-1
    do i=j+1,m
      a(i,j)=(a(i,j)+a(j,i))/2.0E+00
    end do
  end do

  return
end
subroutine fzero ( f, b, c, r, re, ae, iflag )
!
!*******************************************************************************
!
!! FZERO searches for a zero of a function F(X) in a given interval.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    L. f. shampine and h. a. watts,
!    fzero, a root-solving code*,
!    sc-tm-70-631, september 1970.
!
!    t. j. dekker,
!    finding a zero by means of successive linear interpolation*,
!    'constructive aspects of the fundamental theorem of algebra',
!    edited by b. dejon, p. henrici, 1969.
!
!  Discussion:
!
!    fzero searches for a zero of a function f(x) between
!    the given values b and c until the width of the interval
!    (b,c) has collapsed to within a tolerance specified by
!    the stopping criterion, abs(b-c) <= 2.*(rw*abs(b)+ae).
!    the method used is an efficient combination of bisection
!    and the secant rule.
!
!  Parameters:
!
!     f,b,c,r,re and ae are input parameters
!     b,c and iflag are output parameters (flagged by an * below)
!
!        f     - name of the real valued external function.  this name
!                must be in an external statement in the calling
!                program.  f must be a function of one real argument.
!
!       *b     - one end of the interval (b,c).  the value returned for
!                b usually is the better approximation to a zero of f.
!
!       *c     - the other end of the interval (b,c)
!
!        r     - a (better) guess of a zero of f which could help in
!                speeding up convergence.  if f(b) and f(r) have
!                opposite signs, a root will be found in the interval
!                (b,r); if not, but f(r) and f(c) have opposite
!                signs, a root will be found in the interval (r,c);
!                otherwise, the interval (b,c) will be searched for a
!                possible root.  when no better guess is known, it is
!                recommended that r be set to b or c; because if r is
!                not interior to the interval (b,c), it will be ignored.
!
!        re    - relative error used for rw in the stopping criterion.
!                if the requested re is less than machine precision,
!                then rw is set to approximately machine precision.
!
!        ae    - absolute error used in the stopping criterion.  if the
!                given interval (b,c) contains the origin, then a
!                nonzero value should be chosen for ae.
!
!       *iflag - a status code.  user must check iflag after each call.
!                control returns to the user from fzero in all cases.
!
!                1  b is within the requested tolerance of a zero.
!                   the interval (b,c) collapsed to the requested
!                   tolerance, the function changes sign in (b,c), and
!                   f(x) decreased in magnitude as (b,c) collapsed.
!
!                2  f(b) = 0.  however, the interval (b,c) may not have
!                   collapsed to the requested tolerance.
!
!                3  b may be near a singular point of f(x).
!                   the interval (b,c) collapsed to the requested tol-
!                   erance and the function changes sign in (b,c), but
!                   f(x) increased in magnitude as (b,c) collapsed,i.e.
!                     abs(f(b out)) > max(abs(f(b in)),abs(f(c in)))
!
!                4  no change in sign of f(x) was found although the
!                   interval (b,c) collapsed to the requested tolerance.
!                   the user must examine this case and decide whether
!                   b is near a local minimum of f(x), or b is near a
!                   zero of even multiplicity, or neither of these.
!
!                5  too many (> 500) function evaluations used.
!
  real a
  real acbs
  real acmb
  real ae
  real aw
  real b
  real c
  real cmb
  real er
  real f
  real fa
  real fb
  real fc
  real fx
  real fz
  integer ic
  integer iflag
  integer kount
  real p
  real q
  real r
  real re
  real rw
  real t
  real tol
  real z
!
  external f
!
  er = 2.0E+00 * epsilon ( er )
!
!     initialize
!
  z=r
  if ( r<= min (b,c).or.r>= max (b,c)) z=c
  rw= max (re,er)
  aw= max (ae,0.0)
  ic=0
  t=z
  fz=f(t)
  fc=fz
  t=b
  fb=f(t)
  kount=2
  if ( sign(1.0,fz)==sign(1.0,fb)) go to 1
  c=z
  go to 2
    1 if ( z==c) go to 2
  t=c
  fc=f(t)
  kount=3
  if ( sign(1.0,fz)==sign(1.0,fc)) go to 2
  b=z
  fb=fz
    2 a=c
  fa=fc
  acbs=abs(b-c)
  fx= max (abs(fb),abs(fc))
!
    3 if (abs(fc) >= abs(fb)) go to 4
!     perform interchange
  a=b
  fa=fb
  b=c
  fb=fc
  c=a
  fc=fa

4 continue

  cmb=0.5*(c-b)
  acmb=abs(cmb)
  tol=rw*abs(b)+aw
!
!  test stopping criterion and function count
!
  if (acmb <= tol) go to 10
  if ( fb==0.) go to 11
  if ( kount>=500) go to 14
!
!  calculate new iterate implicitly as b+p/q
!     where we arrange p >= 0.
!     the implicit form is used to prevent overflow.
!
  p=(b-a)*fb
  q=fa-fb
  if (p >= 0.) go to 5
  p=-p
  q=-q
!
!  update a and check for satisfactory reduction
!  in the size of the bracketing interval.
!  if not, perform bisection.
!
    5 a=b
  fa=fb
  ic=ic+1
  if (ic < 4) go to 6
  if (8.*acmb >= acbs) go to 8
  ic=0
  acbs=acmb
!
!  test for too small a change
!
    6 if (p > abs(q)*tol) go to 7
!
!  increment by tolerance
!
  b=b+sign(tol,cmb)
  go to 9
!
!  root ought to be between b and (c+b)/2.
!
    7 if (p >= cmb*q) go to 8
!
!  use secant rule
!
  b=b+p/q
  go to 9
!
!  use bisection
!
    8 b=b+cmb
!
!  have completed computation for new iterate b
!
    9 t=b
  fb=f(t)
  kount=kount+1
!
!     decide whether next step is interpolation or extrapolation
!
  if (sign(1.0,fb) /= sign(1.0,fc)) go to 3
  c=a
  fc=fa
  go to 3
!
!
!     finished. process results for proper setting of iflag
!
   10 if (sign(1.0E+00,fb) == sign(1.0E+00,fc)) go to 13
  if (abs(fb) > fx) go to 12
  iflag = 1
  return
   11 iflag = 2
  return
   12 iflag = 3
  return
   13 iflag = 4
  return
   14 iflag = 5
  return
end
subroutine gamlim ( xmin, xmax )
!
!*******************************************************************************
!
!! GAMLIM computes the minimum and maximum bounds for X in GAMMA(X).
!
!
!  Discussion:
!
!    GAMLIM calculates the minimum and maximum legal bounds for X in GAMMA(X).
!
!  Parameters:
!
!    Output, real XMIN, the minimum legal value of X in GAMMA(X).
!    Any smaller value might result in underflow.
!
!    Output, real XMAX, the maximum legal value of X in GAMMA(X).
!    Any larger value will cause overflow.
!
  real alnbig
  real alnsml
  integer i
  real r1mach
  real xln
  real xmax
  real xmin
  real xold
!
  alnsml = log ( tiny ( alnsml ) )
  xmin = -alnsml

  do i=1,10
    xold = xmin
    xln = log(xmin)
    xmin = xmin - xmin*((xmin+0.5)*xln - xmin - 0.2258 + alnsml) &
      / (xmin*xln + 0.5)
    if (abs(xmin-xold)<0.005) go to 20
  end do

  write ( *, * ) ' '
  write ( *, * ) 'GAMLIM - Fatal error!'
  write ( *, * ) '  Unable to determine XMIN.'
  stop

 20   continue

  xmin = -xmin + 0.01

  alnbig = log ( huge ( alnbig ) )
  xmax = alnbig

  do i=1,10
    xold = xmax
    xln = log(xmax)
    xmax = xmax - xmax*((xmax-0.5)*xln - xmax + 0.9189 - alnbig) &
      / (xmax*xln - 0.5)
    if (abs(xmax-xold)<0.005) go to 40
  end do

  write ( *, * ) ' '
  write ( *, * ) 'GAMLIM - Fatal error!'
  write ( *, * ) '  Unable to determine XMAX.'
  stop

 40 continue

  xmax = xmax - 0.01E+00
  xmin = max (xmin, -xmax+1.0E+00 )

  return
end
function gamma ( x )
!
!*******************************************************************************
!
!! GAMMA computes the gamma function.
!
!
!  Parameters:
!
!    Input, real X, the argument of the gamma function, which must not
!    be 0, -1, or any other negative integral value.
!
!    Output, real GAMMA, the value of the gamma function of X.
!
  real csevl
  real dxrel
  real gamma
  real gcs(23)
  integer i
  integer inits
  integer n
  integer ngcs
  real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510
  real r1mach
  real r9lgmc
  real sinpiy
  real sq2pil
  real x
  real xmax
  real xmin
  real y
!
  data gcs( 1) / .008571195590989331e0/
  data gcs( 2) / .004415381324841007e0/
  data gcs( 3) / .05685043681599363e0/
  data gcs( 4) /-.004219835396418561e0/
  data gcs( 5) / .001326808181212460e0/
  data gcs( 6) /-.0001893024529798880e0/
  data gcs( 7) / .0000360692532744124e0/
  data gcs( 8) /-.0000060567619044608e0/
  data gcs( 9) / .0000010558295463022e0/
  data gcs(10) /-.0000001811967365542e0/
  data gcs(11) / .0000000311772496471e0/
  data gcs(12) /-.0000000053542196390e0/
  data gcs(13) / .0000000009193275519e0/
  data gcs(14) /-.0000000001577941280e0/
  data gcs(15) / .0000000000270798062e0/
  data gcs(16) /-.0000000000046468186e0/
  data gcs(17) / .0000000000007973350e0/
  data gcs(18) /-.0000000000001368078e0/
  data gcs(19) / .0000000000000234731e0/
  data gcs(20) /-.0000000000000040274e0/
  data gcs(21) / .0000000000000006910e0/
  data gcs(22) /-.0000000000000001185e0/
  data gcs(23) / .0000000000000000203e0/
!
! sq2pil is log (sqrt (2.*pi) )
!
  data sq2pil / 0.91893853320467274e0/
  data ngcs, xmin, xmax, dxrel /0, 3*0.0E+00 /

!
! initialize.  find legal bounds for x, and determine the number of
! terms in the series required to attain an accuracy ten times better
! than machine precision.
!
  if ( ngcs == 0 ) then
    ngcs = inits (gcs, 23, 0.1*r1mach(3))
    call gamlim (xmin, xmax)
    dxrel = sqrt ( epsilon ( dxrel ) )
  end if

  y = abs(x)
  if (y>10.0) go to 50
!
! compute gamma(x) for abs(x) <= 10.0.  reduce interval and
! find gamma(1+y) for 0. <= y < 1. first of all.
!
  n = x
  if (x<0.0) n = n - 1
  y = x - real(n)
  n = n - 1
  gamma = 0.9375 + csevl(2.*y-1., gcs, ngcs)
  if (n==0) return

  if (n>0) go to 30
!
! compute gamma(x) for x < 1.
!
  n = -n
  if (x==0.) call xerror ( 'gamma   x is 0', 14, 4, 2)

  if (x<0. .and. x+real(n-2)==0.) then
     call xerror (  'gamma   x is a negative integer', 31, 4, 2)
  end if

  if (x<(-0.5) .and. abs((x-aint(x-0.5))/x)<dxrel) then
    call xerror ( &
      'gamma   answer lt half precision because x too near negative integer', &
      68, 1, 1)
  end if

  do i=1,n
    gamma = gamma / (x+real(i-1))
  end do

  return
!
! gamma(x) for x >= 2.
!
 30 continue

  do i=1,n
    gamma = (y+real(i))*gamma
  end do

  return
!
! compute gamma(x) for abs(x) > 10.0.  recall y = abs(x).
!
 50   continue

  if (x>xmax) then
    call xerror ( 'gamma   x so big gamma overflows', 32, 3, 2)
  end if

  gamma = 0.
  if (x<xmin) then
    call xerror ( 'gamma   x so small gamma underflows', 35, 2, 1)
  end if

  if (x<xmin) return

  gamma = exp((y-0.5)*log(y) - y + sq2pil + r9lgmc(y) )
  if (x>0.) return

  if (abs((x-aint(x-0.5))/x)<dxrel) then
    call xerror ( &
      'gamma   answer lt half precision, x too near negative integer', &
      61, 1, 1)
  end if

  sinpiy = sin (pi*y)

  if (sinpiy==0.) then
    write ( *, * ) ' '
    write ( *, * ) 'GAMMA - Fatal error!'
    write ( *, * ) '  The argument X is a negative integer.'
    stop
  end if

  gamma = -pi / (y*sinpiy*gamma)

  return
end
subroutine gl15t ( f, a, b, xl, xr, r, ae, ra, rasc, fmin, fmax )
!
!*******************************************************************************
!
!! GL15T computes integral of g(x) over (a,b), with error estimate.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!            on entry
!              f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs
!                       to be declared e x t e r n a l in the
!                       calling program.
!                       the function g(x) is defined to be
!                       g(x)=f(phi(x))*phip(x)
!                       where phi(x) is the cubic given by
!                       the arithmetic statement function below.
!                       phip(x) is its derivative.  the variables
!                       xl and xr are the left and right endpoints
!                       of a parent interval of which (a,b) is a part.
!
!              a      - real
!                       lower limit of integration
!
!              b      - real
!                       upper limit of integration
!
!              xl     - double precision
!              xr     - double precision
!                       lower and upper limits of parent interval
!                       of which [a,b] is a part.
!
!            on return
!              r - real
!                       approximation to the integral i
!                       r is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal
!                       addition of abscissae to the 7-point gauss
!                       rule (resg).
!
!              ae - real
!                       estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-r)
!
!              ra - real
!                       approximation to the integral j
!
!              rasc - real
!                       approximation to the integral of abs(g-i/(b-a))
!                       over (a,b)
!
!              fmax, fmin - real
!                       max and min values of the function f on (a,b)
!
  real a
  double precision absc
  real ae
  real b
  double precision centr
  real dhlgth
  real, save :: epmach = 0.0E+00
  real f
  real fc
  real fmax
  real fmin
  real fsum
  real fv1(7)
  real fv2(7)
  real fval1
  real fval2
  real hlgth
  integer j
  integer jtw
  integer jtwm1
  real phi
  real phip
  real phiu
  real r
  real r1mach
  real ra
  real rasc
  real resg
  real resk
  real reskh
  real sl
  real sr
  double precision u
  real, save :: uflow = 0.0E+00
  real wg(4)
  real wgk(8)
  real xgk(8)
  double precision xl
  double precision xr
!
  external f
!
!           the abscissae and weights are given for the interval (-1,1)
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!           wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule
!
  data wg(  1) / 0.129484966168869693270611432679082d0 /
  data wg(  2) / 0.279705391489276667901467771423780d0 /
  data wg(  3) / 0.381830050505118944950369775488975d0 /
  data wg(  4) / 0.417959183673469387755102040816327d0 /
!
  data xgk(  1) / 0.991455371120812639206854697526329d0 /
  data xgk(  2) / 0.949107912342758524526189684047851d0 /
  data xgk(  3) / 0.864864423359769072789712788640926d0 /
  data xgk(  4) / 0.741531185599394439863864773280788d0 /
  data xgk(  5) / 0.586087235467691130294144838258730d0 /
  data xgk(  6) / 0.405845151377397166906606412076961d0 /
  data xgk(  7) / 0.207784955007898467600689403773245d0 /
  data xgk(  8) / 0.000000000000000000000000000000000d0 /
!
  data wgk(  1) / 0.022935322010529224963732008058970d0 /
  data wgk(  2) / 0.063092092629978553290700663189204d0 /
  data wgk(  3) / 0.104790010322250183839876322541518d0 /
  data wgk(  4) / 0.140653259715525918745189590510238d0 /
  data wgk(  5) / 0.169004726639267902826583426598550d0 /
  data wgk(  6) / 0.190350578064785409913256402421014d0 /
  data wgk(  7) / 0.204432940075298892414161999234649d0 /
  data wgk(  8) / 0.209482141084727828012999174891714d0 /
!
  phi(u)=xr-(xr-xl)*u*u*(2.*u+3.)
  phip(u)=-6.*u*(u+1.)
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - r of the 7-point gauss formula
!           resk   - r of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  if ( epmach == 0.0E+00 ) then
    epmach = epsilon ( 1.0E+00 )
    uflow = tiny ( uflow )
  end if
!
  if ( xl<xr ) then
     sl=sngl(xl)
     sr=sngl(xr)
  else
     sl=sngl(xr)
     sr=sngl(xl)
  end if

  hlgth = 0.5 * (b-a)
  centr = a+hlgth
  dhlgth = abs(hlgth)
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the absolute error.
!
  u=(centr-xr)/(xr-xl)
  phiu=phi(u)
  if ( phiu<=sl .or. phiu>=sr) phiu=centr
  fmin=f(phiu)
  fmax=fmin
  fc=fmin*phip(u)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  ra = abs(resk)

  do j=1,3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    u=(centr-absc-xr)/(xr-xl)
    phiu=phi(u)
    if ( phiu<=sl .or. phiu>=sr) phiu=centr
    fval1=f(phiu)
    fmax=max(fmax,fval1)
    fmin=min(fmin,fval1)
    fval1=fval1*phip(u)
    u=(centr+absc-xr)/(xr-xl)
    phiu=phi(u)
    if ( phiu<=sl .or. phiu>=sr) phiu=centr
    fval2=f(phiu)
    fmax=max(fmax,fval2)
    fmin=min(fmin,fval2)
    fval2=fval2*phip(u)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    ra = ra+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1,4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    u=(centr-absc-xr)/(xr-xl)
    phiu=phi(u)
    if ( phiu<=sl .or. phiu>=sr) phiu=centr
    fval1=f(phiu)
    fmax=max(fmax,fval1)
    fmin=min(fmin,fval1)
    fval1=fval1*phip(u)
    u=(centr+absc-xr)/(xr-xl)
    phiu=phi(u)
    if ( phiu<=sl .or. phiu>=sr) phiu=centr
    fval2=f(phiu)
    fmax=max(fmax,fval2)
    fmin=min(fmin,fval2)
    fval2=fval2*phip(u)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    ra = ra+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk*0.5
  rasc = wgk(8)*abs(fc-reskh)

  do j=1,7
    rasc = rasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  r = resk*hlgth
  ra = ra*dhlgth
  rasc = rasc*dhlgth
  ae = abs((resk-resg)*hlgth)

  if ( rasc/=0.0.and.ae/=0.0) then
    ae = rasc*min(1.0,(0.2e+03*ae/rasc)**1.5)
  end if

  if ( ra>uflow/(0.5e+02*epmach)) then
    ae = max((epmach*0.5e+02)*ra,ae)
  end if

  return
end
subroutine grdchk ( n, x, fcn, f, g, typsiz, sx, fscale, rnf, analtl, &
  wrk1, msg, ipr )
!
!*******************************************************************************
!
!! GRDCHK checks an analytic gradient against an estimated gradient.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
! n            --> dimension of problem
! x(n)         --> estimate to a root of fcn
! fcn          --> name of subroutine to evaluate optimization function
!                  must be declared external in calling routine
!                       fcn:  r(n) --> r(1)
! f            --> function value:  fcn(x)
! g(n)         --> gradient:  g(x)
! typsiz(n)    --> typical size for each component of x
! sx(n)        --> diagonal scaling matrix:  sx(i)=1./typsiz(i)
! fscale       --> estimate of scale of objective function fcn
! rnf          --> relative noise in optimization function fcn
! analtl       --> tolerance for comparison of estimated and
!                  analytical gradients
! wrk1(n)      --> workspace
! msg         <--  message or error code
!                    on output: =-21, probable coding error of gradient
! ipr          --> device to which to send output
!
  integer n
!
  real analtl
  real f
  real fscale
  real g(n)
  real gs
  integer i
  integer ipr
  integer ker
  integer msg
  real rnf
  real sx(n)
  real typsiz(n)
  real wrk
  real wrk1(n)
  real x(n)
!
  external fcn
!
!  Compute the first order finite difference gradient;
!  compare it to the analytic gradient.
!
  call fstofd(1,1,n,x,fcn,f,wrk1,sx,rnf,wrk,1)

  ker = 0

  do i = 1, n

    gs = max ( abs(f), fscale ) / max ( abs(x(i)), typsiz(i) )

    if ( abs(g(i)-wrk1(i))>max(abs(g(i)),gs)*analtl) then
      ker=1
    end if

  end do

  if ( ker /= 0 ) then
    write ( ipr, * ) ' '
    write ( ipr, * ) 'GRDCHK - probable error in analytic gradient.'
    write ( ipr, * ) ' '
    write ( ipr, * ) ' grdchk     comp            analytic            est'
    write ( ipr, 902 ) (i,g(i),wrk1(i),i=1,n)
    msg = -21
  end if

  return

  902 format(' grdchk    ',i5,3x,e20.13,3x,e20.13)
end
subroutine heschk ( nr, n, x, fcn, d1fcn, d2fcn, f, g, a, typsiz, sx, rnf, &
  analtl, iagflg, udiag, wrk1, wrk2, msg, ipr )
!
!*******************************************************************************
!
!! HESCHK checks an analytic hessian against a computed estimate.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> estimate to a root of fcn
! fcn          --> name of subroutine to evaluate optimization function
!                  must be declared external in calling routine
!                       fcn:  r(n) --> r(1)
! d1fcn        --> name of subroutine to evaluate gradient of fcn.
!                  must be declared external in calling routine
! d2fcn        --> name of subroutine to evaluate hessian of fcn.
!                  must be declared external in calling routine
! f            --> function value:  fcn(x)
! g(n)        <--  gradient:  g(x)
! a(n,n)      <--  on exit:  hessian in lower triangular part and diag
! typsiz(n)    --> typical size for each component of x
! sx(n)        --> diagonal scaling matrix:  sx(i)=1./typsiz(i)
! rnf          --> relative noise in optimization function fcn
! analtl       --> tolerance for comparison of estimated and
!                  analytical gradients
! iagflg       --> =1 if analytic gradient supplied
! udiag(n)     --> workspace
! wrk1(n)      --> workspace
! wrk2(n)      --> workspace
! msg         <--> message or error code
!                    on input : if =1xx do not compare anal + est hess
!                    on output: =-22, probable coding error of hessian
! ipr          --> device to which to send output
!
  integer n
  integer nr
!
  real a(nr,n)
  real analtl
  real f
  real g(n)
  real hs
  integer i
  integer iagflg
  integer ipr
  integer j
  integer ker
  integer msg
  real rnf
  real sx(n)
  real typsiz(n)
  real udiag(n)
  real wrk1(n)
  real wrk2(n)
  real x(n)
!
  external fcn
  external d1fcn
  external d2fcn
!
!  Compute finite difference approximation a to the hessian.
!
  if ( iagflg==1) then
    call fstofd(nr,n,n,x,d1fcn,g,a,sx,rnf,wrk1,3)
  else
    call sndofd(nr,n,x,fcn,f,a,sx,rnf,wrk1,wrk2)
  end if

  ker=0
!
! copy lower triangular part of "a" to upper triangular part
! and diagonal of "a" to udiag
!
  do j=1,n
    udiag(j)=a(j,j)
    do i=j+1,n
      a(j,i)=a(i,j)
    end do
  end do
!
! compute analytic hessian and compare to finite difference
! approximation.
!
  call d2fcn(nr,n,x,a)

  do j=1,n

    hs=max(abs(g(j)),1.0)/max(abs(x(j)),typsiz(j))

    if ( abs(a(j,j)-udiag(j))>max(abs(udiag(j)),hs)*analtl) then
      ker=1
    end if

    do i=j+1,n
      if ( abs(a(i,j)-a(j,i))>max(abs(a(i,j)),hs)*analtl) then
        ker=1
      end if
    end do

  end do

  if ( ker==0) go to 90

    write(ipr,901)
    do i=1,n
      do j=1,i-1
        write(ipr,902) i,j,a(i,j),a(j,i)
      end do
      write(ipr,902) i,i,a(i,i),udiag(i)
    end do

    msg=-22
!     end if
   90 continue
  return
  901 format('heschk    probable error in coding of analytic hessian.'/ &
             'heschk      row  col',14x,'analytic',14x,'(estimate)')
  902 format('heschk    ',2i5,2x,e20.13,2x,'(',e20.13,')')
end
subroutine hookdr ( nr, n, x, f, g, a, udiag, p, xpls, fpls, fcn, sx, stepmx, &
  steptl, dlt, iretcd, mxtake, amu, dltp, phi, phip0, sc, xplsp, wrk0, epsm, &
  itncnt, ipr )
!
!*******************************************************************************
!
!! HOOKDR finds the next Newton iterate by the More-Hebdon method.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> old iterate x[k-1]
! f            --> function value at old iterate, f(x)
! g(n)         --> gradient at old iterate, g(x), or approximate
! a(n,n)       --> cholesky decomposition of hessian in lower
!                  triangular part and diagonal.
!                  hessian in upper triangular part and udiag.
! udiag(n)     --> diagonal of hessian in a(.,.)
! p(n)         --> newton step
! xpls(n)     <--  new iterate x[k]
! fpls        <--function value at new iterate, f(xpls)
! fcn          --> name of subroutine to evaluate function
! sx(n)        --> diagonal scaling matrix for x
! stepmx       --> maximum allowable step size
! steptl       --> relative step size at which successive iterates
!                  considered close enough to terminate algorithm
! dlt         <--> trust region radius
! iretcd      <--  return code
!                    =0 satisfactory xpls found
!                    =1 failed to find satisfactory xpls sufficiently
!                       distinct from x
! mxtake      <--  boolean flag indicating step of maximum length used
! amu         <--> [retain value between successive calls]
! dltp        <--> [retain value between successive calls]
! phi         <--> [retain value between successive calls]
! phip0       <--> [retain value between successive calls]
! sc(n)        --> workspace
! xplsp(n)     --> workspace
! wrk0(n)      --> workspace
! epsm         --> machine epsilon
! itncnt       --> iteration count
! ipr          --> device to which to send output
!
  integer n
  integer nr
!
  real a(nr,n)
  real alpha
  real amu
  real beta
  real dlt
  real dltp
  real epsm
  real f
  real fpls
  real fplsp
  logical fstime
  real g(n)
  integer i
  integer ipr
  integer iretcd
  integer itncnt
  integer j
  logical mxtake
  logical nwtake
  real p(n)
  real phi
  real phip0
  real rnwtln
  real sc(n)
  real stepmx
  real steptl
  real sx(n)
  real tmp
  real udiag(n)
  real x(n)
  real xpls(n)
  real xplsp(n)
  real wrk0(n)
!
  external fcn
!
  iretcd=4
  fstime=.true.

  tmp=0.
  do i=1,n
    tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
  end do

  rnwtln=sqrt(tmp)

  if ( itncnt>1) go to 100
    amu=0.
!
!  if first iteration and trust region not provided by user,
!  compute initial trust region.
!
    if ( dlt/= (-1.)) go to 100

      alpha=0.
      do i=1,n
        alpha=alpha+(g(i)*g(i))/(sx(i)*sx(i))
      end do

      beta=0.0E+00
      do i=1,n
        tmp=0.0E+00
        do j=i,n
          tmp=tmp + (a(j,i)*g(j))/(sx(j)*sx(j))
        end do
        beta=beta+tmp*tmp
      end do
      dlt=alpha*sqrt(alpha)/beta
      dlt = min(dlt, stepmx)
!       end if
!     end if
!
  100 continue
!
!  Find new step by more-hebdon algorithm.
!
  call hookst(nr,n,g,a,udiag,p,sx,rnwtln,dlt,amu,dltp,phi,phip0, &
    fstime,sc,nwtake,wrk0,epsm,ipr)

  dltp=dlt
!
!  check new point and update trust region
!
  call tregup(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl, &
    dlt,iretcd,xplsp,fplsp,xpls,fpls,mxtake,ipr,3,udiag)

  if ( iretcd > 1 ) then
    go to 100
  end if

  return
end
subroutine hookst ( nr, n, g, a, udiag, p, sx, rnwtln, dlt, amu, &
  dltp, phi, phip0, fstime, sc, nwtake, wrk0, epsm, ipr )
!
!*******************************************************************************
!
!! HOOKST finds the new step by the More-Hebdon algorithm.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! g(n)         --> gradient at current iterate, g(x)
! a(n,n)       --> cholesky decomposition of hessian in
!                  lower triangular part and diagonal.
!                  hessian or approx in upper triangular part
! udiag(n)     --> diagonal of hessian in a(.,.)
! p(n)         --> newton step
! sx(n)        --> diagonal scaling matrix for n
! rnwtln       --> newton step length
! dlt         <--> trust region radius
! amu         <--> [retain value between successive calls]
! dltp         --> trust region radius at last exit from this routine
! phi         <--> [retain value between successive calls]
! phip0       <--> [retain value between successive calls]
! fstime      <--> boolean. =.true. if first entry to this routine
!                  during k-th iteration
! sc(n)       <--  current step
! nwtake      <--  boolean, =.true. if newton step taken
! wrk0         --> workspace
! epsm         --> machine epsilon
! ipr          --> device to which to send output
!
  integer n
  integer nr
!
  real a(nr,n)
  real addmax
  real alo
  real amu
  real amulo
  real amuup
  real dlt
  real dltp
  logical done
  real epsm
  logical fstime
  real g(n)
  real hi
  integer i
  integer ipr
  integer j
  logical nwtake
  real p(n)
  real phi
  real phip
  real phip0
  real rnwtln
  real sc(n)
  real snrm2
  real stepln
  real sx(n)
  real udiag(n)
  real wrk0(n)
!
! hi and alo are constants used in this routine.
! change here if other values are to be substituted.
!
  ipr=ipr
  hi=1.5
  alo=.75
  if ( rnwtln>hi*dlt) go to 15
!     if ( rnwtln<=hi*dlt)
!     then
!
!       take newton step
!
    nwtake=.true.
    do i=1,n
      sc(i)=p(i)
    end do
    dlt=min(dlt,rnwtln)
    amu=0.0E+00
    return

!     else
!
!       newton step not taken
!
   15   continue
    nwtake=.false.
    if ( amu<=0.) go to 20
!       if ( amu>0.)
!       then
      amu=amu- (phi+dltp) *((dltp-dlt)+phi)/(dlt*phip)
!       end if
   20   continue
    phi=rnwtln-dlt
    if ( .not.fstime) go to 28
!       if ( fstime)
!       then
      wrk0(1:n) = sx(1:n) * sx(1:n) * p(1:n)
!
!  solve l*y = (sx**2)*p
!
      call forslv(nr,n,a,wrk0,wrk0)

      phip0=-snrm2(n,wrk0,1)**2/rnwtln
      fstime=.false.

   28   phip = phip0
    amulo = -phi/phip
    amuup = 0.0E+00
    do i = 1, n
      amuup = amuup+(g(i)*g(i))/(sx(i)*sx(i))
    end do
    amuup = sqrt(amuup)/dlt
    done=.false.
!
!  test value of amu; generate next amu if necessary
!
  100   continue
    if ( done) return
    if ( amu>=amulo .and. amu<=amuup) go to 110
!       if ( amu<amulo .or.  amu>amuup)
!       then
      amu=max(sqrt(amulo*amuup),amuup*1.0e-3)
!       end if
  110   continue
!
!  copy (h,udiag) to l
!  where h <-- h+amu*(sx**2) [do not actually change (h,udiag)]
!
    do j=1,n

      a(j,j)=udiag(j) + amu*sx(j)*sx(j)
      a(j+1:n,j)=a(j,j+1:n)

    end do
!
!  factor h=l(l+)
!
    call choldc(nr,n,a,0.0,sqrt(epsm),addmax)
!
!  solve h*p = l(l+)*sc = -g
!
    wrk0(1:n)=-g(1:n)

    call lltslv(nr,n,a,sc,wrk0)
!
!  reset h.  note since udiag has not been destroyed we need do
!  nothing here.  h is in the upper part and in udiag, still intact
!
    stepln=0.0E+00
    do i=1,n
      stepln=stepln + sx(i)*sx(i)*sc(i)*sc(i)
    end do

    stepln=sqrt(stepln)
    phi=stepln-dlt

    do i=1,n
      wrk0(i)=sx(i)*sx(i)*sc(i)
    end do

    call forslv(nr,n,a,wrk0,wrk0)

    phip=-snrm2(n,wrk0,1)**2/stepln

    if ( (alo*dlt>stepln .or. stepln>hi*dlt) .and.(amuup-amulo>0.)) then
      go to 170
    end if

!       if ( (alo*dlt<=stepln .and. stepln<=hi*dlt) .or.
!            (amuup-amulo<=0.))
!       then
!
!         sc is acceptable hookstep
!
      done=.true.
      go to 100
!       else
!
!         sc not acceptable hookstep.  select new amu
!
  170     continue
      amulo=max(amulo,amu-(phi/phip))
      if ( phi<0.) amuup=min(amuup,amu)
      amu=amu-(stepln*phi)/(dlt*phip)
      go to 100
!       end if
!     end if
!
  951 format('0hookst    take newton step')
  952 format('0hookst    newton step not taken')
  953 format(' hookst    sc is not acceptable')
  954 format(' hookst    sc is acceptable')
  955 format(' hookst    current step (sc)')
  956 format(' hookst    amu   =',e20.13)
  957 format(' hookst    amulo =',e20.13)
  958 format(' hookst    amuup =',e20.13)
  959 format(' hookst    phi   =',e20.13)
  960 format(' hookst    phip  =',e20.13)
  961 format(' hookst    dlt   =',e20.13/' hookst    stepln=',e20.13)
  962 format('0hookst    find new amu')
  963 format(' hookst       ',5(e20.13,3x))
end
subroutine hsnint ( nr, n, a, sx, method )
!
!*******************************************************************************
!
!! HSNINT provides initial hessian when using secant updates.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(n,n)      <--  initial hessian (lower triangular matrix)
! sx(n)        --> diagonal scaling matrix for x
! method       --> algorithm to use to solve minimization problem
!                    =1,2 factored secant method used
!                    =3   unfactored secant method used
!
  integer n
  integer nr
!
  real a(nr,n)
  integer j
  integer method
  real sx(n)
!
  do j = 1, n

    if ( method == 3 ) then
      a(j,j)=sx(j)**2
    else
      a(j,j)=sx(j)
    end if

    a(j+1:n,j)=0.0E+00

  end do

  return
end
function i1mach ( i )
!
!*******************************************************************************
!
!! I1MACH returns integer machine constants.
!
!
!  I/O unit numbers.
!
!    I1MACH(1) = the standard input unit.
!    I1MACH(2) = the standard output unit.
!    I1MACH(3) = the standard punch unit.
!    I1MACH(4) = the standard error message unit.
!
!  Words.
!
!    I1MACH(5) = the number of bits per integer storage unit.
!    I1MACH(6) = the number of characters per integer storage unit.
!
!  Integers.
!
!  Assume integers are represented in the S digit base A form:
!
!  Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
!  where 0<=X(I)<A for I=0 to S-1.
!
!    I1MACH(7) = A, the base.
!    I1MACH(8) = S, the number of base A digits.
!    I1MACH(9) = A**S-1, the largest integer.
!
!  Floating point numbers
!
!  Assume floating point numbers are represented in the T digit base B form:
!
!    Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
!
!  where 0<=X(I)<B for I=1 to T, 0<X(1) and EMIN<=E<=EMAX
!
!    I1MACH(10) = B, the base.
!
!  Single precision
!
!    I1MACH(11) = T, the number of base B digits.
!    I1MACH(12) = EMIN, the smallest exponent E.
!    I1MACH(13) = EMAX, the largest exponent E.
!
!  Double precision
!
!    I1MACH(14) = T, the number of base B digits.
!    I1MACH(15) = EMIN, the smallest exponent E.
!    I1MACH(16) = EMAX, the largest exponent E.
!
!  To alter this function for a particular environment, the desired set of DATA
!  statements should be activated by removing the C from column 1.  On rare
!  machines, a STATIC statement may need to be added, but probably more systems
!  prohibit than require it.
!
!  Also, the values of I1MACH(1) through I1MACH(4) should be checked for
!  consistency with the local operating system.  For FORTRAN 77, you may wish
!  to adjust the data statement so imach(6) is set to 1, and then to comment
!  out the executable test on I.EQ.6 below.
!
!  For IEEE-arithmetic machines (binary standard), the first set of constants
!  below should be appropriate, except perhaps for IMACH(1) - IMACH(4).
!
  integer i
  integer i1mach
  integer imach(16)
  integer output
!
  equivalence (imach(4),output)
!
!  IEEE arithmetic machines, such as the ATT 3B series, Motorola
!  68000 based machines such as the SUN 3 and ATT PC 7300, and
!  8087 based micros such asthe IBM PC and ATT 6300.
!
   data imach( 1) /    5 /
   data imach( 2) /    6 /
   data imach( 3) /    7 /
   data imach( 4) /    6 /
   data imach( 5) /   32 /
   data imach( 6) /    4 /
   data imach( 7) /    2 /
   data imach( 8) /   31 /
   data imach( 9) / 2147483647 /
   data imach(10) /    2 /
   data imach(11) /   24 /
   data imach(12) / -125 /
   data imach(13) /  128 /
   data imach(14) /   53 /
   data imach(15) / -1021 /
   data imach(16) /  1024 /
!
!  ALLIANT FX/8 UNIX FORTRAN compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  AMDAHL machines.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  BURROUGHS 1700 system.
!
!      data imach( 1) /    7 /
!      data imach( 2) /    2 /
!      data imach( 3) /    2 /
!      data imach( 4) /    2 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   33 /
!      data imach( 9) / Z1FFFFFFFF /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -256 /
!      data imach(13) /  255 /
!      data imach(14) /   60 /
!      data imach(15) / -256 /
!      data imach(16) /  255 /
!
!  BURROUGHS 5700 system.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -50 /
!      data imach(16) /  76 /
!
!  BURROUGHS 6700/7700 systems.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  48 /
!      data imach( 6) /   6 /
!      data imach( 7) /   2 /
!      data imach( 8) /  39 /
!      data imach( 9) / O0007777777777777 /
!      data imach(10) /   8 /
!      data imach(11) /  13 /
!      data imach(12) / -50 /
!      data imach(13) /  76 /
!      data imach(14) /  26 /
!      data imach(15) / -32754 /
!      data imach(16) /  32780 /
!
!  CDC CYBER 170/180 series using NOS
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   48 /
!      data imach(12) / -974 /
!      data imach(13) / 1070 /
!      data imach(14) /   96 /
!      data imach(15) / -927 /
!      data imach(16) / 1070 /
!
!  CDC CYBER 170/180 series using NOS/VE
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     7 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) / 9223372036854775807 /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -4095 /
!      data imach(13) /  4094 /
!      data imach(14) /    94 /
!      data imach(15) / -4095 /
!      data imach(16) /  4094 /
!
!  CDC CYBER 200 series
!
!      data imach( 1) /      5 /
!      data imach( 2) /      6 /
!      data imach( 3) /      7 /
!      data imach( 4) /      6 /
!      data imach( 5) /     64 /
!      data imach( 6) /      8 /
!      data imach( 7) /      2 /
!      data imach( 8) /     47 /
!      data imach( 9) / X'00007FFFFFFFFFFF' /
!      data imach(10) /      2 /
!      data imach(11) /     47 /
!      data imach(12) / -28625 /
!      data imach(13) /  28718 /
!      data imach(14) /     94 /
!      data imach(15) / -28625 /
!      data imach(16) /  28718 /
!
!  CDC 6000/7000 series using FTN4.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / 00007777777777777777B /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CDC 6000/7000 series using FTN5.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   60 /
!      data imach( 6) /   10 /
!      data imach( 7) /    2 /
!      data imach( 8) /   48 /
!      data imach( 9) / O"00007777777777777777" /
!      data imach(10) /    2 /
!      data imach(11) /   47 /
!      data imach(12) / -929 /
!      data imach(13) / 1070 /
!      data imach(14) /   94 /
!      data imach(15) / -929 /
!      data imach(16) / 1069 /
!
!  CONVEX C-1.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  CONVEX C-120 (native mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (native mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1023 /
!      data imach(13) /  1023 /
!      data imach(14) /    53 /
!      data imach(15) / -1023 /
!      data imach(16) /  1023 /
!
!  CONVEX C-120 (IEEE mode) without -R8 option
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    0 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CONVEX C-120 (IEEE mode) with -R8 option
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     0 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    53 /
!      data imach(12) / -1021 /
!      data imach(13) /  1024 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  CRAY 1, 2, XMP and YMP.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /   102 /
!      data imach( 4) /     6 /
!      data imach( 5) /    64 /
!      data imach( 6) /     8 /
!      data imach( 7) /     2 /
!      data imach( 8) /    63 /
!      data imach( 9) /  777777777777777777777B /
!      data imach(10) /     2 /
!      data imach(11) /    47 /
!      data imach(12) / -8189 /
!      data imach(13) /  8190 /
!      data imach(14) /    94 /
!      data imach(15) / -8099 /
!      data imach(16) /  8190 /
!
!  DATA GENERAL ECLIPSE S/200.
!
!      data imach( 1) /   11 /
!      data imach( 2) /   12 /
!      data imach( 3) /    8 /
!      data imach( 4) /   10 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) /32767 /
!      data imach(10) /   16 /
!      data imach(11) /    6 /
!      data imach(12) /  -64 /
!      data imach(13) /   63 /
!      data imach(14) /   14 /
!      data imach(15) /  -64 /
!      data imach(16) /   63 /
!
!  ELXSI 6400
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     6 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  HARRIS 220
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HARRIS SLASH 6 and SLASH 7.
!
!      data imach( 1) /       5 /
!      data imach( 2) /       6 /
!      data imach( 3) /       0 /
!      data imach( 4) /       6 /
!      data imach( 5) /      24 /
!      data imach( 6) /       3 /
!      data imach( 7) /       2 /
!      data imach( 8) /      23 /
!      data imach( 9) / 8388607 /
!      data imach(10) /       2 /
!      data imach(11) /      23 /
!      data imach(12) /    -127 /
!      data imach(13) /     127 /
!      data imach(14) /      38 /
!      data imach(15) /    -127 /
!      data imach(16) /     127 /
!
!  HONEYWELL DPS 8/70 and 600/6000 series.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /   43 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   63 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  HP 2100, 3 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   39 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 2100, 4 word double precision option with FTN4
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    4 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   23 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   55 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  HP 9000
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     7 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1015 /
!      data imach(16) /  1017 /
!
!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!  and PERKIN ELMER (INTERDATA) 3230.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   7 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z7FFFFFFF /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  63 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  63 /
!
!  IBM PC - Microsoft FORTRAN
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!
!      data imach( 1) /     4 /
!      data imach( 2) /     7 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   127 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1023 /
!
!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!  For the INTERDATA FORTRAN VII compiler, replace the Z's specifying hex
!  constants with Y's.
!
!      data imach( 1) /   5 /
!      data imach( 2) /   6 /
!      data imach( 3) /   6 /
!      data imach( 4) /   6 /
!      data imach( 5) /  32 /
!      data imach( 6) /   4 /
!      data imach( 7) /   2 /
!      data imach( 8) /  31 /
!      data imach( 9) / Z'7FFFFFFF' /
!      data imach(10) /  16 /
!      data imach(11) /   6 /
!      data imach(12) / -64 /
!      data imach(13) /  62 /
!      data imach(14) /  14 /
!      data imach(15) / -64 /
!      data imach(16) /  62 /
!
!  PDP-10 (KA processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   54 /
!      data imach(15) / -101 /
!      data imach(16) /  127 /
!
!  PDP-10 (KI processor).
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    5 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / "377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   62 /
!      data imach(15) / -128 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 32-bit integer arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PDP-11 FORTRANS supporting 16-bit integer arithmetic.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  PRIME 50 series systems with 32-bit integers and 64V MODE instructions,
!  supplied by Igor Bray.
!
!      data imach( 1) /            1 /
!      data imach( 2) /            1 /
!      data imach( 3) /            2 /
!      data imach( 4) /            1 /
!      data imach( 5) /           32 /
!      data imach( 6) /            4 /
!      data imach( 7) /            2 /
!      data imach( 8) /           31 /
!      data imach( 9) / :17777777777 /
!      data imach(10) /            2 /
!      data imach(11) /           23 /
!      data imach(12) /         -127 /
!      data imach(13) /         +127 /
!      data imach(14) /           47 /
!      data imach(15) /       -32895 /
!      data imach(16) /       +32637 /
!
!  SEQUENT BALANCE 8000.
!
!      data imach( 1) /     0 /
!      data imach( 2) /     0 /
!      data imach( 3) /     7 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     1 /
!      data imach( 7) /     2 /
!      data imach( 8) /    31 /
!      data imach( 9) /  2147483647 /
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -125 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  SUN Microsystems UNIX F77 compiler.
!
!      data imach( 1) /     5 /
!      data imach( 2) /     6 /
!      data imach( 3) /     6 /
!      data imach( 4) /     0 /
!      data imach( 5) /    32 /
!      data imach( 6) /     4 /
!      data imach( 7) /     2 /
!      data imach( 8) /    32 /
!      data imach( 9) /2147483647/
!      data imach(10) /     2 /
!      data imach(11) /    24 /
!      data imach(12) /  -126 /
!      data imach(13) /   128 /
!      data imach(14) /    53 /
!      data imach(15) / -1022 /
!      data imach(16) /  1024 /
!
!  SUN 3 (68881 or FPA)
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    6 /
!      data imach( 4) /    0 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -125 /
!      data imach(13) /  128 /
!      data imach(14) /   53 /
!      data imach(15) / -1021 /
!      data imach(16) /  1024 /
!
!  UNIVAC 1100 series.
!  Note that the punch unit, I1MACH(3), has been set to 7, which is appropriate
!  for the UNIVAC-FOR system.  If you have the UNIVAC-FTN system, set it to 1
!  instead.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   36 /
!      data imach( 6) /    6 /
!      data imach( 7) /    2 /
!      data imach( 8) /   35 /
!      data imach( 9) / O377777777777 /
!      data imach(10) /    2 /
!      data imach(11) /   27 /
!      data imach(12) / -128 /
!      data imach(13) /  127 /
!      data imach(14) /   60 /
!      data imach(15) /-1024 /
!      data imach(16) / 1023 /
!
!  VAX.
!
!      data imach( 1) /    5 /
!      data imach( 2) /    6 /
!      data imach( 3) /    7 /
!      data imach( 4) /    6 /
!      data imach( 5) /   32 /
!      data imach( 6) /    4 /
!      data imach( 7) /    2 /
!      data imach( 8) /   31 /
!      data imach( 9) / 2147483647 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
!  Z80 microprocessor.
!
!      data imach( 1) /    1 /
!      data imach( 2) /    1 /
!      data imach( 3) /    0 /
!      data imach( 4) /    1 /
!      data imach( 5) /   16 /
!      data imach( 6) /    2 /
!      data imach( 7) /    2 /
!      data imach( 8) /   15 /
!      data imach( 9) / 32767 /
!      data imach(10) /    2 /
!      data imach(11) /   24 /
!      data imach(12) / -127 /
!      data imach(13) /  127 /
!      data imach(14) /   56 /
!      data imach(15) / -127 /
!      data imach(16) /  127 /
!
  if ( i < 1 .or. i > 16 )then
    write ( *, * ) ' '
    write(*,*)'I1MACH - Fatal error!'
    write(*,*)'I is out of bounds:',i
    i1mach=0
    stop
  else
    i1mach=imach(i)
  end if

  return
end
subroutine i_factor ( n, ifac )
!
!*******************************************************************************
!
!! I_FACTOR factors an integer.
!
!
!  Modified:
!
!    14 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the number to be factored.
!
!    Output, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer i
  integer ib
  integer ifac(15)
  integer j
  integer n
  integer nf
  integer nl
  integer nq
  integer nr
  integer ntry
!
  ifac(1) = n

  nf = 0
  nl = n

  if ( n == 0 ) then
    nf = 1
    ifac(2) = nf
    ifac(2+nf) = 0
    return
  end if

  if ( n < 1 ) then
    nf = nf + 1
    ifac(2+nf) = -1
    nl = - n
  end if

  if ( nl == 1 ) then
    nf = nf + 1
    ifac(2) = nf
    ifac(2+nf) = 1
    return
  end if

  j = 0

  do while ( nl > 1 )

    j = j + 1
!
!  Choose a trial divisor, NTRY.
!
    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if
!
!  Divide by the divisor as many times as possible.
!
    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nl = nq
      nf = nf + 1
!
!  Make sure factors of 2 appear in the front of the list.
!
      if ( ntry /= 2 ) then

        ifac(2+nf) = ntry

      else

        do i = nf, 2, -1
          ifac(i+2) = ifac(i+1)
        end do
        ifac(3) = 2

      end if

    end do

  end do

  ifac(2) = nf

  return
end
function inits ( os, nos, eta )
!
!*******************************************************************************
!
!! INITS estimates the order of an orthogonal series guaranteeing a given accuracy.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Ordinarily, ETA will be chosen to be one-tenth machine precision.
!
!             input arguments --
! os     array of nos coefficients in an orthogonal series.
! nos    number of coefficients in os.
! eta    requested accuracy of series.
!
  integer nos

  integer ii,i
  integer inits
  real os,err,eta
  dimension os(nos)
!
  if (nos<1) then
    call xerror ( 'inits   number of coefficients lt 1', 35, 2, 2)
  end if

  err = 0.
  do ii=1,nos
    i = nos + 1 - ii
    err = err + abs(os(i))
    if (err>eta) go to 20
  end do

 20   if (i==nos) call xerror ( 'inits   eta may be too small', 28, 1, 2)
  inits = i

  return
end
function isamax ( n, x, incx )
!
!*******************************************************************************
!
!! ISAMAX finds the index of the vector element of maximum absolute value.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of SX.
!
!    Output, integer ISAMAX, the index of the element of SX of maximum
!    absolute value.
!
  integer i
  integer incx
  integer isamax
  integer ix
  integer n
  real samax
  real x(*)
!
  if ( n <= 0 ) then

    isamax = 0

  else if ( n == 1 ) then

    isamax = 1

  else if ( incx == 1 ) then

    isamax = 1
    samax = abs ( x(1) )

    do i = 2, n

      if ( abs ( x(i) ) > samax ) then
        isamax = i
        samax = abs ( x(i) )
      end if

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    isamax = 1
    samax = abs ( x(ix) )

    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > samax ) then
        isamax = i
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function j4save ( iwhich, ivalue, iset )
!
!*******************************************************************************
!
!! J4SAVE saves variables needed by the library error handling routines.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters
!      --input--
!        iwhich - index of item desired.
!                = 1 refers to current error number.
!                = 2 refers to current error control flag.
!                 = 3 refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                 = 4 refers to the maximum number of times any
!                     message is to be printed (as set by xermax).
!                 = 5 refers to the total number of units to which
!                     each error message is to be written.
!                 = 6 refers to the 2nd unit for error messages
!                 = 7 refers to the 3rd unit for error messages
!                 = 8 refers to the 4th unit for error messages
!                 = 9 refers to the 5th unit for error messages
!        ivalue - the value to be set for the iwhich-th parameter,
!                 if iset is .true. .
!        iset   - if iset=.true., the iwhich-th parameter will be
!                 given the value, ivalue.  if iset=.false., the
!                 iwhich-th parameter will be unchanged, and ivalue
!                 is a dummy parameter.
!      --output--
!        the (old) value of the iwhich-th parameter will be returned
!        in the function value, j4save.
!
!     written by ron jones, with slatec common math library subcommittee
!    adapted from bell laboratories port library error handler
!     latest revision: 23 may 1979
!
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer, save, dimension ( 9 ) :: iparam = (/ 0, 2, 0, 10, 1, 0, 0, 0, 0 /)
  logical iset
  integer ivalue
  integer iwhich
  integer j4save
!
  j4save = iparam(iwhich)
  if (iset) then
    iparam(iwhich) = ivalue
  end if

  return
end
function inbin ( x, nbins, xmin, xmax, width )
!
!*******************************************************************************
!
!! INBIN takes a real value X and finds the correct bin for it.
!
!
!  Discussion:
!
!    Values below XMIN come back in 1.  Values above XMAX come back
!    in NBINS.
!
!  Parameters:
!
!    Input, real X, a value to be binned.
!
!    Input, integer NBINS, the number of bins.
!
!    Input, real XMIN, XMAX, the minimum and maximum bin limits.
!
!    Input, real WIDTH, the width of each bin.
!
!    Output, integer INBIN, the index of the bin containing X.
!
  integer inbin
  integer nbins
  real width
  real x
  real xmax
  real xmin
!
  if ( x < xmin ) then
    inbin = 1
  else if ( x >= xmax ) then
    inbin = nbins
  else
    inbin = 2 + int ( ( x - xmin ) / width )
  end if

  return
end
subroutine jairy ( x, rx, c, ai, dai )
!
!*******************************************************************************
!
!! JAIRY computes the Airy function and its derivative.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!                                   input
!
!         x - argument, computed by asyjy, x unrestricted
!        rx - rx=sqrt(abs(x)), computed by asyjy
!         c - c=2.*(abs(x)**1.5)/3., computed by asyjy
!
!                                  output
!
!        ai - value of function ai(x)
!       dai - value of the derivative dai(x)
!
!                                written by
!
!                                d. e. amos
!                               s. l. daniel
!                               m. k. weston
!
  real a
  real ai
  real ajn
  real ajp, ak1, ak2, ak3
  real b(15)
  real c, ccv, con2, con3
  real con4, con5, cv
  real da
  real dai, dajn, dajp, dak1, dak2, dak3, db, ec
  real e1, e2, fpi12, f1, f2
  integer i
  integer j
  integer m1
  integer m1d
  integer m2
  integer m2d
  integer m3
  integer m3d
  integer m4
  integer m4d
  integer n1
  integer n1d
  integer n2
  integer n2d, n3, n3d, n4, n4d
  real rtrx, rx, scv, t, temp1, temp2, tt
  real x
  dimension ajp(19), ajn(19), a(15)
  dimension ak1(14), ak2(23), ak3(14)
  dimension dajp(19), dajn(19), da(15), db(15)
  dimension dak1(14), dak2(24), dak3(14)

  data n1,n2,n3,n4/14,23,19,15/
  data m1,m2,m3,m4/12,21,17,13/
  data fpi12,con2,con3,con4,con5/ &
    1.30899693899575, 5.03154716196777, &
    3.80004589867293e-01, 8.33333333333333e-01, 8.66025403784439e-01/

  data ak1(1), ak1(2), ak1(3), ak1(4), ak1(5), ak1(6), ak1(7), &
          ak1(8), ak1(9), ak1(10),ak1(11),ak1(12),ak1(13), &
          ak1(14)         / 2.20423090987793e-01,-1.25290242787700e-01, &
      1.03881163359194e-02, 8.22844152006343e-04,-2.34614345891226e-04, &
      1.63824280172116e-05, 3.06902589573189e-07,-1.29621999359332e-07, &
      8.22908158823668e-09, 1.53963968623298e-11,-3.39165465615682e-11, &
      2.03253257423626e-12,-1.10679546097884e-14,-5.16169497785080e-15/

  data ak2(1), ak2(2), ak2(3), ak2(4), ak2(5), ak2(6), ak2(7), &
          ak2(8), ak2(9), ak2(10),ak2(11),ak2(12),ak2(13),ak2(14), &
          ak2(15),ak2(16),ak2(17),ak2(18),ak2(19),ak2(20),ak2(21), &
          ak2(22),ak2(23) / 2.74366150869598e-01, 5.39790969736903e-03, &
     -1.57339220621190e-03, 4.27427528248750e-04,-1.12124917399925e-04, &
      2.88763171318904e-05,-7.36804225370554e-06, 1.87290209741024e-06, &
     -4.75892793962291e-07, 1.21130416955909e-07,-3.09245374270614e-08, &
      7.92454705282654e-09,-2.03902447167914e-09, 5.26863056595742e-10, &
     -1.36704767639569e-10, 3.56141039013708e-11,-9.31388296548430e-12, &
      2.44464450473635e-12,-6.43840261990955e-13, 1.70106030559349e-13, &
     -4.50760104503281e-14, 1.19774799164811e-14,-3.19077040865066e-15/

  data ak3(1), ak3(2), ak3(3), ak3(4), ak3(5), ak3(6), ak3(7), &
          ak3(8), ak3(9), ak3(10),ak3(11),ak3(12),ak3(13), &
          ak3(14)         / 2.80271447340791e-01,-1.78127042844379e-03, &
      4.03422579628999e-05,-1.63249965269003e-06, 9.21181482476768e-08, &
     -6.52294330229155e-09, 5.47138404576546e-10,-5.24408251800260e-11, &
      5.60477904117209e-12,-6.56375244639313e-13, 8.31285761966247e-14, &
     -1.12705134691063e-14, 1.62267976598129e-15,-2.46480324312426e-16/

  data ajp(1), ajp(2), ajp(3), ajp(4), ajp(5), ajp(6), ajp(7), &
          ajp(8), ajp(9), ajp(10),ajp(11),ajp(12),ajp(13),ajp(14), &
          ajp(15),ajp(16),ajp(17),ajp(18), &
          ajp(19)         / 7.78952966437581e-02,-1.84356363456801e-01, &
      3.01412605216174e-02, 3.05342724277608e-02,-4.95424702513079e-03, &
     -1.72749552563952e-03, 2.43137637839190e-04, 5.04564777517082e-05, &
     -6.16316582695208e-06,-9.03986745510768e-07, 9.70243778355884e-08, &
      1.09639453305205e-08,-1.04716330588766e-09,-9.60359441344646e-11, &
      8.25358789454134e-12, 6.36123439018768e-13,-4.96629614116015e-14, &
     -3.29810288929615e-15, 2.35798252031104e-16/

  data ajn(1), ajn(2), ajn(3), ajn(4), ajn(5), ajn(6), ajn(7), &
          ajn(8), ajn(9), ajn(10),ajn(11),ajn(12),ajn(13),ajn(14), &
          ajn(15),ajn(16),ajn(17),ajn(18), &
          ajn(19)         / 3.80497887617242e-02,-2.45319541845546e-01, &
      1.65820623702696e-01, 7.49330045818789e-02,-2.63476288106641e-02, &
     -5.92535597304981e-03, 1.44744409589804e-03, 2.18311831322215e-04, &
     -4.10662077680304e-05,-4.66874994171766e-06, 7.15218807277160e-07, &
      6.52964770854633e-08,-8.44284027565946e-09,-6.44186158976978e-10, &
      7.20802286505285e-11, 4.72465431717846e-12,-4.66022632547045e-13, &
     -2.67762710389189e-14, 2.36161316570019e-15/

  data a(1),   a(2),   a(3),   a(4),   a(5),   a(6),   a(7), &
          a(8),   a(9),   a(10),  a(11),  a(12),  a(13),  a(14), &
          a(15)           / 4.90275424742791e-01, 1.57647277946204e-03, &
     -9.66195963140306e-05, 1.35916080268815e-07, 2.98157342654859e-07, &
     -1.86824767559979e-08,-1.03685737667141e-09, 3.28660818434328e-10, &
     -2.57091410632780e-11,-2.32357655300677e-12, 9.57523279048255e-13, &
     -1.20340828049719e-13,-2.90907716770715e-15, 4.55656454580149e-15, &
     -9.99003874810259e-16/

  data b(1),   b(2),   b(3),   b(4),   b(5),   b(6),   b(7), &
          b(8),   b(9),   b(10),  b(11),  b(12),  b(13),  b(14), &
          b(15)           / 2.78593552803079e-01,-3.52915691882584e-03, &
     -2.31149677384994e-05, 4.71317842263560e-06,-1.12415907931333e-07, &
     -2.00100301184339e-08, 2.60948075302193e-09,-3.55098136101216e-11, &
     -3.50849978423875e-11, 5.83007187954202e-12,-2.04644828753326e-13, &
     -1.10529179476742e-13, 2.87724778038775e-14,-2.88205111009939e-15, &
     -3.32656311696166e-16/

  data n1d,n2d,n3d,n4d/14,24,19,15/
  data m1d,m2d,m3d,m4d/12,22,17,13/

  data dak1(1), dak1(2), dak1(3), dak1(4), dak1(5), dak1(6), &
          dak1(7), dak1(8), dak1(9), dak1(10),dak1(11),dak1(12), &
         dak1(13),dak1(14)/ 2.04567842307887e-01,-6.61322739905664e-02, &
     -8.49845800989287e-03, 3.12183491556289e-03,-2.70016489829432e-04, &
     -6.35636298679387e-06, 3.02397712409509e-06,-2.18311195330088e-07, &
     -5.36194289332826e-10, 1.13098035622310e-09,-7.43023834629073e-11, &
      4.28804170826891e-13, 2.23810925754539e-13,-1.39140135641182e-14/

  data dak2(1), dak2(2), dak2(3), dak2(4), dak2(5), dak2(6), &
          dak2(7), dak2(8), dak2(9), dak2(10),dak2(11),dak2(12), &
          dak2(13),dak2(14),dak2(15),dak2(16),dak2(17),dak2(18), &
          dak2(19),dak2(20),dak2(21),dak2(22),dak2(23), &
          dak2(24)        / 2.93332343883230e-01,-8.06196784743112e-03, &
      2.42540172333140e-03,-6.82297548850235e-04, 1.85786427751181e-04, &
     -4.97457447684059e-05, 1.32090681239497e-05,-3.49528240444943e-06, &
      9.24362451078835e-07,-2.44732671521867e-07, 6.49307837648910e-08, &
     -1.72717621501538e-08, 4.60725763604656e-09,-1.23249055291550e-09, &
      3.30620409488102e-10,-8.89252099772401e-11, 2.39773319878298e-11, &
     -6.48013921153450e-12, 1.75510132023731e-12,-4.76303829833637e-13, &
      1.29498241100810e-13,-3.52679622210430e-14, 9.62005151585923e-15, &
     -2.62786914342292e-15/

  data dak3(1), dak3(2), dak3(3), dak3(4), dak3(5), dak3(6), &
          dak3(7), dak3(8), dak3(9), dak3(10),dak3(11),dak3(12), &
         dak3(13),dak3(14)/ 2.84675828811349e-01, 2.53073072619080e-03, &
     -4.83481130337976e-05, 1.84907283946343e-06,-1.01418491178576e-07, &
      7.05925634457153e-09,-5.85325291400382e-10, 5.56357688831339e-11, &
     -5.90889094779500e-12, 6.88574353784436e-13,-8.68588256452194e-14, &
      1.17374762617213e-14,-1.68523146510923e-15, 2.55374773097056e-16/

  data dajp(1), dajp(2), dajp(3), dajp(4), dajp(5), dajp(6), &
          dajp(7), dajp(8), dajp(9), dajp(10),dajp(11),dajp(12), &
          dajp(13),dajp(14),dajp(15),dajp(16),dajp(17),dajp(18), &
          dajp(19)        / 6.53219131311457e-02,-1.20262933688823e-01, &
      9.78010236263823e-03, 1.67948429230505e-02,-1.97146140182132e-03, &
     -8.45560295098867e-04, 9.42889620701976e-05, 2.25827860945475e-05, &
     -2.29067870915987e-06,-3.76343991136919e-07, 3.45663933559565e-08, &
      4.29611332003007e-09,-3.58673691214989e-10,-3.57245881361895e-11, &
      2.72696091066336e-12, 2.26120653095771e-13,-1.58763205238303e-14, &
     -1.12604374485125e-15, 7.31327529515367e-17/

  data dajn(1), dajn(2), dajn(3), dajn(4), dajn(5), dajn(6), &
          dajn(7), dajn(8), dajn(9), dajn(10),dajn(11),dajn(12), &
          dajn(13),dajn(14),dajn(15),dajn(16),dajn(17),dajn(18), &
          dajn(19)        / 1.08594539632967e-02, 8.53313194857091e-02, &
     -3.15277068113058e-01,-8.78420725294257e-02, 5.53251906976048e-02, &
      9.41674060503241e-03,-3.32187026018996e-03,-4.11157343156826e-04, &
      1.01297326891346e-04, 9.87633682208396e-06,-1.87312969812393e-06, &
     -1.50798500131468e-07, 2.32687669525394e-08, 1.59599917419225e-09, &
     -2.07665922668385e-10,-1.24103350500302e-11, 1.39631765331043e-12, &
      7.39400971155740e-14,-7.32887475627500e-15/

  data da(1),  da(2),  da(3),  da(4),  da(5),  da(6),  da(7), &
          da(8),  da(9),  da(10), da(11), da(12), da(13), da(14), &
          da(15)          / 4.91627321104601e-01, 3.11164930427489e-03, &
      8.23140762854081e-05,-4.61769776172142e-06,-6.13158880534626e-08, &
      2.87295804656520e-08,-1.81959715372117e-09,-1.44752826642035e-10, &
      4.53724043420422e-11,-3.99655065847223e-12,-3.24089119830323e-13, &
      1.62098952568741e-13,-2.40765247974057e-14, 1.69384811284491e-16, &
      8.17900786477396e-16/

  data db(1),  db(2),  db(3),  db(4),  db(5),  db(6),  db(7), &
          db(8),  db(9),  db(10), db(11), db(12), db(13), db(14), &
          db(15)          /-2.77571356944231e-01, 4.44212833419920e-03, &
     -8.42328522190089e-05,-2.58040318418710e-06, 3.42389720217621e-07, &
     -6.24286894709776e-09,-2.36377836844577e-09, 3.16991042656673e-10, &
     -4.40995691658191e-12,-5.18674221093575e-12, 9.64874015137022e-13, &
     -4.90190576608710e-14,-1.77253430678112e-14, 5.55950610442662e-15, &
     -7.11793337579530e-16/
!
  if (x<0.0) go to 90
  if (c>5.0) go to 60
  if (x>1.20) go to 30
  t = (x+x-1.2)*con4
  tt = t + t
  j = n1
  f1 = ak1(j)
  f2 = 0.0E+00

  do i=1,m1
    j = j - 1
    temp1 = f1
    f1 = tt*f1 - f2 + ak1(j)
    f2 = temp1
  end do

  ai = t*f1 - f2 + ak1(1)

  j = n1d
  f1 = dak1(j)
  f2 = 0.0E+00
  do i=1,m1d
    j = j - 1
    temp1 = f1
    f1 = tt*f1 - f2 + dak1(j)
    f2 = temp1
  end do

  dai = -(t*f1-f2+dak1(1))

  return
!
   30 continue
  t = (x+x-con2)*con3
  tt = t + t
  j = n2
  f1 = ak2(j)
  f2 = 0.0E+00

  do i=1,m2
    j = j - 1
    temp1 = f1
    f1 = tt*f1 - f2 + ak2(j)
    f2 = temp1
  end do

  rtrx = sqrt(rx)
  ec = exp(-c)
  ai = ec*(t*f1-f2+ak2(1))/rtrx
  j = n2d
  f1 = dak2(j)
  f2 = 0.0E+00

  do i =1,m2d
    j = j - 1
    temp1 = f1
    f1 = tt*f1 - f2 + dak2(j)
    f2 = temp1
  end do
  dai = -ec*(t*f1-f2+dak2(1))*rtrx
  return

   60 continue
  t = 10.0/c - 1.0E+00
  tt = t + t
  j = n1
  f1 = ak3(j)
  f2 = 0.0E+00

  do i=1,m1
    j = j - 1
    temp1 = f1
    f1 = tt*f1 - f2 + ak3(j)
    f2 = temp1
  end do
  rtrx = sqrt(rx)
  ec = exp(-c)
  ai = ec*(t*f1-f2+ak3(1))/rtrx
  j = n1d
  f1 = dak3(j)
  f2 = 0.0E+00

  do i=1,m1d
    j = j - 1
    temp1 = f1
    f1 = tt*f1 - f2 + dak3(j)
    f2 = temp1
  end do
  dai = -rtrx*ec*(t*f1-f2+dak3(1))
  return

   90 continue
  if (c>5.0) go to 120
  t = 0.4*c - 1.0E+00
  tt = t + t
  j = n3
  f1 = ajp(j)
  e1 = ajn(j)
  f2 = 0.0E+00
  e2 = 0.0E+00

  do i=1,m3
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt*f1 - f2 + ajp(j)
    e1 = tt*e1 - e2 + ajn(j)
    f2 = temp1
    e2 = temp2
  end do

  ai = (t*e1-e2+ajn(1)) - x*(t*f1-f2+ajp(1))
  j = n3d
  f1 = dajp(j)
  e1 = dajn(j)
  f2 = 0.0E+00
  e2 = 0.0E+00

  do i=1,m3d
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt*f1 - f2 + dajp(j)
    e1 = tt*e1 - e2 + dajn(j)
    f2 = temp1
    e2 = temp2
  end do

  dai = x*x*(t*f1-f2+dajp(1)) + (t*e1-e2+dajn(1))

  return
!
  120 continue
  t = 10.0/c - 1.0E+00
  tt = t + t
  j = n4
  f1 = a(j)
  e1 = b(j)
  f2 = 0.0E+00
  e2 = 0.0E+00

  do i=1,m4
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt*f1 - f2 + a(j)
    e1 = tt*e1 - e2 + b(j)
    f2 = temp1
    e2 = temp2
  end do
  temp1 = t*f1 - f2 + a(1)
  temp2 = t*e1 - e2 + b(1)
  rtrx = sqrt(rx)
  cv = c - fpi12
  ccv = cos(cv)
  scv = sin(cv)
  ai = (temp1*ccv-temp2*scv)/rtrx
  j = n4d
  f1 = da(j)
  e1 = db(j)
  f2 = 0.0E+00
  e2 = 0.0E+00

  do i = 1, m4d
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt*f1 - f2 + da(j)
    e1 = tt*e1 - e2 + db(j)
    f2 = temp1
    e2 = temp2
  end do

  temp1 = t*f1 - f2 + da(1)
  temp2 = t*e1 - e2 + db(1)
  e1 = ccv*con5 + 0.5*scv
  e2 = scv*con5 - 0.5*ccv
  dai = (temp1*e1-temp2*e2)*rtrx

  return
end
subroutine lltslv ( nr, n, a, x, b )
!
!*******************************************************************************
!
!! LLTSLV solves A*x=b where A = L * Transpose(L), and only L is stored.
!
!
!  Discussion:
!
!    L is a lower triangular matrix.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(n,n)       --> matrix of form l(l-transpose).
!                  on return a is unchanged.
! x(n)        <--  solution vector
! b(n)         --> right-hand side vector
!                  if b is not required by calling program, then
!                  b and x may share the same storage.
!
  integer n
  integer nr
!
  real a(nr,n)
  real b(n)
  real x(n)
!
!  Forward solve, result in X.
!
  call forslv ( nr, n, a, x, b )
!
!  Back solve, result in X.
!
  call bakslv ( nr, n, a, x, x )

  return
end
subroutine lnsrch ( n, x, f, g, p, xpls, fpls, fcn, mxtake, iretcd, stepmx, &
  steptl, sx, ipr )
!
!*******************************************************************************
!
!! LNSRCH finds a next newton iterate by line search.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! n            --> dimension of problem
! x(n)         --> old iterate:   x[k-1]
! f            --> function value at old iterate, f(x)
! g(n)         --> gradient at old iterate, g(x), or approximate
! p(n)         --> non-zero newton step
! xpls(n)     <--  new iterate x[k]
! fpls        <--function value at new iterate, f(xpls)
! fcn          --> name of subroutine to evaluate function
! iretcd      <--  return code
! mxtake      <--  boolean flag indicating step of maximum length used
! stepmx       --> maximum allowable step size
! steptl       --> relative step size at which successive iterates
!                  considered close enough to terminate algorithm
! sx(n)        --> diagonal scaling matrix for x
! ipr          --> device to which to send output
!
! internal variables
!
! sln              newton length
! rln              relative length of newton step
!
  integer n
!
  real a
  real almbda
  real b
  real disc
  real f
  real fpls
  real g(n)
  integer i
  integer ipr
  integer iretcd
  logical mxtake
  real p(n)
  real pfpls
  real plmbda
  real rln
  real rmnlmb
  real scl
  real sln
  real slp
  real stepmx
  real steptl
  real sx(n)
  real t1
  real t2
  real t3
  real tlmbda
  real tmp
  real x(n)
  real xpls(n)
!
  external fcn
!
  mxtake = .false.
  iretcd = 2

  sln = sqrt ( sum ( ( sx(1:n) * p(1:n) )**2 ) )
!
!  Newton step longer than maximum allowed.
!
  if ( sln > stepmx ) then
    scl=stepmx / sln
    p(1:n) = scl * p(1:n)
    sln = stepmx
  end if

  slp = dot_product ( g, p )

  rln=0.
  do i=1,n
    rln=max(rln,abs(p(i))/max(abs(x(i)),1./sx(i)))
  end do

  rmnlmb=steptl/rln
  almbda=1.0E+00
!
! loop
! check if new iterate satisfactory.  generate new lambda if necessary.
!
  100 continue

  if ( iretcd < 2 ) return

  xpls(1:n) = x(1:n) + almbda * p(1:n)

  call fcn ( n, xpls, fpls )

  if ( fpls> f+slp*1.0e-4*almbda) go to 130
!
! solution found
!
    iretcd=0
    if ( almbda==1.0E+00 .and. sln> .99*stepmx) mxtake=.true.
    go to 100
!
! solution not (yet) found
!
!     else
  130   if ( almbda >= rmnlmb) go to 140
!       if ( almbda < rmnlmb)
!       then
!
! no satisfactory xpls found sufficiently distinct from x
!
      iretcd=1
      go to 100
!
! calculate new lambda
!
  140     continue

          if ( almbda/=1.0) go to 150
!         if ( almbda==1.0)
!         then
!
! first backtrack: quadratic fit
!
        tlmbda=-slp/(2.*(fpls-f-slp))
        go to 170
!         else
!
! all subsequent backtracks: cubic fit
!
  150       t1=fpls-f-almbda*slp
        t2=pfpls-f-plmbda*slp
        t3=1.0/(almbda-plmbda)
        a=t3*(t1/(almbda*almbda) - t2/(plmbda*plmbda))
        b=t3*(t2*almbda/(plmbda*plmbda) - t1*plmbda/(almbda*almbda) )
        disc=b*b-3.0*a*slp
        if ( disc<= b*b) go to 160
!           if ( disc> b*b)
!           then
!
! only one positive critical point, must be minimum
!
          tlmbda=(-b+sign(1.0,a)*sqrt(disc))/(3.0*a)
          go to 165
!           else
!
! both critical points positive, first is minimum
!
  160         tlmbda=(-b-sign(1.0,a)*sqrt(disc))/(3.0*a)
!           end if
  165       if ( tlmbda> .5*almbda) tlmbda=.5*almbda
!         end if
  170     plmbda=almbda
      pfpls=fpls
      if ( tlmbda>= almbda*.1) go to 180
!         if ( tlmbda<almbda/10.)
!         then
        almbda=almbda*.1
        go to 190
!         else
  180       almbda=tlmbda
!         end if
!       end if
!     end if
  190 go to 100
  950 format(' lnsrch    almbda=',e20.13)
  951 format(' lnsrch    new iterate (xpls)')
  952 format(' lnsrch    sln   =',e20.13/' lnsrch    slp   =',e20.13/ &
             ' lnsrch    rmnlmb=',e20.13/' lnsrch    stepmx=',e20.13/ &
             ' lnsrch    steptl=',e20.13)
  953 format(' lnsrch    f(xpls)=',e20.13)
  954 format('0lnsrch    newton step (p)')
  955 format(' lnsrch       ',5(e20.13,3x))
end
subroutine mvmltl ( nr, n, a, x, y )
!
!*******************************************************************************
!
!! MVMLTL computes y = L*x where l is a lower triangular matrix stored in A.
!
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(n,n)       --> lower triangular (n*n) matrix
! x(n)         --> operand vector
! y(n)        <--  result vector
!
! note
!
! x and y cannot share storage
!
  integer n
  integer nr
!
  real a(nr,n)
  integer i
  integer j
  real sum
  real x(n)
  real y(n)
!
  do i = 1, n
    sum = 0.0E+00
    do j = 1, i
      sum = sum + a(i,j) * x(j)
    end do
    y(i) = sum
  end do

  return
end
subroutine mvmlts ( nr, n, a, x, y )
!
!*******************************************************************************
!
!! MVMLTS computes y = A*x where "a" is a symmetric (n*n) matrix stored in its lower
! triangular part and x,y are n-vectors
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(n,n)       --> symmetric (n*n) matrix stored in
!                  lower triangular part and diagonal
! x(n)         --> operand vector
! y(n)        <--  result vector
!
! note
!
! x and y cannot share storage.
!
  integer n
  integer nr
!
  real a(nr,n)
  integer i
  integer j
  real sum
  real x(n)
  real y(n)
!
  do i = 1, n

    sum = 0.0E+00

    do j = 1, i
      sum = sum + a(i,j) * x(j)
    end do

    do j = i+1, n
      sum = sum + a(j,i) * x(j)
    end do

    y(i) = sum

  end do

  return
end
subroutine mvmltu ( nr, n, a, x, y )
!
!*******************************************************************************
!
!! MVMLTU computes y=(l+)x where l is a lower triangular matrix stored in a
! (l-transpose (l+) is taken implicitly)
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
! parameters
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! a(nr,1)       --> lower triangular (n*n) matrix
! x(n)         --> operand vector
! y(n)        <--  result vector
!
! note
!
! x and y cannot share storage
!
  integer n
  integer nr
!
  real a(nr,1)
  integer i
  integer j
  real sum
  real x(n)
  real y(n)
!
  do i = 1, n
    sum = 0.0E+00
    do j = i, n
      sum = sum + a(j,i) * x(j)
    end do
    y(i) = sum
  end do

  return
end
function numxer ( nerr )
!
!*******************************************************************************
!
!! NUMXER returns the most recent error number.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision:  7 june 1978
!
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
  integer j4save
  integer nerr
  integer numxer
!
  nerr = j4save(1,0,.false.)
  numxer = nerr
  return
end
subroutine optchk ( n, x, typsiz, sx, fscale, gradtl, itnlim, ndigit, epsm, &
  dlt, method, iexp, iagflg, iahflg, stepmx, msg, ipr )
!
!*******************************************************************************
!
!! OPTCHK checks the input to the optimization routine.
!
!
!  Parameters:
!
! n            --> dimension of problem
! x(n)         --> on entry, estimate to root of fcn
! typsiz(n)   <--> typical size of each component of x
! sx(n)       <--  diagonal scaling matrix for x
! fscale      <--> estimate of scale of objective function fcn
! gradtl       --> tolerance at which gradient considered close
!                  enough to zero to terminate algorithm
! itnlim      <--> maximum number of allowable iterations
! ndigit      <--> number of good digits in optimization function fcn
! epsm         --> machine epsilon
! dlt         <--> trust region radius
! method      <--> algorithm indicator
! iexp        <--> expense flag
! iagflg      <--> =1 if analytic gradient supplied
! iahflg      <--> =1 if analytic hessian supplied
! stepmx      <--> maximum step size
! msg         <--> message and error code
! ipr          --> device to which to send output
!
  integer n
!
  real dlt
  real epsm
  real fscale
  real gradtl
  integer i
  integer iagflg
  integer iahflg
  integer iexp
  integer ipr
  integer itnlim
  integer method
  integer msg
  integer ndigit
  real stepmx
  real stpsiz
  real sx(n)
  real typsiz(n)
  real x(n)
!
! check that parameters only take on acceptable values.
! if not, set them to default values.
!
  if ( method<1 .or. method>3) method=1
  if ( iagflg/=1) iagflg=0
  if ( iahflg/=1) iahflg=0
  if ( iexp/=0) iexp=1
  if ( mod(msg/2,2)==1 .and. iagflg==0) go to 830
  if ( mod(msg/4,2)==1 .and. iahflg==0) go to 835
!
!  check dimension of problem
!
  if ( n<=0) go to 805
  if ( n==1 .and. mod(msg,2)==0) go to 810
!
!  compute scale matrix
!
  do i=1,n
    if ( typsiz(i)==0.) typsiz(i)=1.0E+00
    if ( typsiz(i)<0.) typsiz(i)=-typsiz(i)
    sx(i)=1.0/typsiz(i)
  end do
!
!  check maximum step size
!
  if (stepmx > 0.0) go to 20
  stpsiz = 0.0E+00
  do i = 1, n
     stpsiz = stpsiz + x(i)*x(i)*sx(i)*sx(i)
  end do

  stpsiz = sqrt(stpsiz)

  stepmx = max(1.0e3*stpsiz, 1.0e3)
   20 continue
!
! check function scale
!
  if ( fscale==0.) fscale=1.0E+00
  if ( fscale<0.) fscale=-fscale
!
! check gradient tolerance
!
  if ( gradtl<0.) go to 815
!
! check iteration limit
!
  if ( itnlim<=0) go to 820
!
! check number of digits of accuracy in function fcn
!
  if ( ndigit==0) go to 825
  if ( ndigit<0) ndigit=-log10(epsm)
!
! check trust region radius
!
  if ( dlt<=0.) dlt=-1.0E+00
  if (dlt > stepmx) dlt = stepmx
  return
!
! error exits
!
  805 write(ipr,901) n
  msg=-1
  go to 895
  810 write(ipr,902)
  msg=-2
  go to 895
  815 write(ipr,903) gradtl
  msg=-3
  go to 895
  820 write(ipr,904) itnlim
  msg=-4
  go to 895
  825 write(ipr,905) ndigit
  msg=-5
  go to 895
  830 write(ipr,906) msg,iagflg
  msg=-6
  go to 895
  835 write(ipr,907) msg,iahflg
  msg=-7
  895 return
  901 format('0optchk    illegal dimension, n=',i5)
  902 format('0optchk    +++ warning +++  this package is inefficient', &
    'for problems of size n=1.'/ &
    ' optchk    check installation libraries for more appropriate routines.'/ &
    ' optchk    if none, set msg and resubmit.')
  903 format('0optchk    illegal tolerance.  gradtl=',e20.13)
  904 format('0optchk    illegal iteration limit.  itnlim=',i5)
  905 format('0optchk    minimization function has no good digits.', &
     'ndigit=',i5)
  906 format('0optchk    user requests that analytic gradient be', &
     ' accepted as properly coded (msg=',i5, '),'/ &
     ' optchk    but analytic gradient not supplied (iagflg=',i5, ').')
  907 format('0optchk    user requests that analytic hessian be', &
     ' accepted as properly coded (msg=',i5, '),'/ &
     ' optchk    but analytic hessian not supplied(iahflg=',i5, ').')
end
subroutine optdrv(nr,n,x,fcn,d1fcn,d2fcn,typsiz,fscale,method,iexp,msg, &
  ndigit,itnlim,iagflg,iahflg,ipr,dlt,gradtl,stepmx,steptl,xpls,fpls,gpls, &
  itrmcd,a,udiag,g,p,sx,wrk0,wrk1,wrk2,wrk3)
!
!*******************************************************************************
!
!! OPTDRV is a driver for the nonlinear optimization package.
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> on entry: estimate to a root of fcn
! fcn          --> name of subroutine to evaluate optimization function
!                  must be declared external in calling routine
!                            fcn: r(n) --> r(1)
! d1fcn        --> (optional) name of subroutine to evaluate gradient
!                  of fcn.  must be declared external in calling routine
! d2fcn        --> (optional) name of subroutine to evaluate hessian of
!                  of fcn.  must be declared external in calling routine
! typsiz(n)    --> typical size for each component of x
! fscale       --> estimate of scale of objective function
! method       --> algorithm to use to solve minimization problem
!                    =1 line search
!                    =2 double dogleg
!                    =3 more-hebdon
! iexp         --> =1 if optimization function fcn is expensive to
!                  evaluate, =0 otherwise.  if set then hessian will
!                  be evaluated by secant update instead of
!                  analytically or by finite differences
! msg         <--> on input:  (>0) message to inhibit certain
!                    automatic checks
!                  on output: (<0) error code; =0 no error
! ndigit       --> number of good digits in optimization function fcn
! itnlim       --> maximum number of allowable iterations
! iagflg       --> =1 if analytic gradient supplied
! iahflg       --> =1 if analytic hessian supplied
! ipr          --> device to which to send output
! dlt          --> trust region radius
! gradtl       --> tolerance at which gradient considered close
!                  enough to zero to terminate algorithm
! stepmx       --> maximum allowable step size
! steptl       --> relative step size at which successive iterates
!                  considered close enough to terminate algorithm
! xpls(n)     <--> on exit:  xpls is local minimum
! fpls        <--> on exit:function value at solution, xpls
! gpls(n)     <--> on exit:  gradient at solution xpls
! itrmcd      <--  termination code
! a(n,n)       --> workspace for hessian (or estimate)
!                  and its cholesky decomposition
! udiag(n)     --> workspace [for diagonal of hessian]
! g(n)         --> workspace (for gradient at current iterate)
! p(n)         --> workspace for step
! sx(n)        --> workspace (for diagonal scaling matrix)
! wrk0(n)      --> workspace
! wrk1(n)      --> workspace
! wrk2(n)      --> workspace
! wrk3(n)      --> workspace
!
!
! internal variables
!
! analtl           tolerance for comparison of estimated and
!                  analytical gradients and hessians
! epsm             machine epsilon
! f              function value: fcn(x)
! itncnt           current iteration, k
! rnf              relative noise in optimization function fcn.
!                       noise=10.**(-ndigit)
!
  integer n
  integer nr
!
  real a(nr,n)
  real amu
  real amusav
  real analtl
  real dlpsav
  real dlt
  real dltp
  real dltsav
  real epsm
  real f
  real fpls
  real fscale
  real g(n)
  real gpls(n)
  real gradtl
  integer i
  integer iagflg
  integer iahflg
  integer icscmx
  integer iexp
  integer ipr
  integer iretcd
  integer itncnt
  integer itnlim
  integer itrmcd
  integer method
  integer msg
  logical mxtake
  integer ndigit
  logical noupdt
  real p(n)
  real phi
  real phip0
  real phisav
  real phpsav
  real r1mach
  real rnf
  real stepmx
  real steptl
  real sx(n)
  real typsiz(n)
  real udiag(n)
  real wrk
  real wrk0(n)
  real wrk1(n)
  real wrk2(n)
  real wrk3(n)
  real x(n)
  real xpls(n)
!
  external fcn,d1fcn,d2fcn
!
! initialization
!
  p(1:n) = 0.0E+00
  itncnt=0
  iretcd=-1
  epsm=r1mach(4)

  call optchk(n,x,typsiz,sx,fscale,gradtl,itnlim,ndigit,epsm, &
    dlt,method,iexp,iagflg,iahflg,stepmx,msg,ipr)
  if ( msg<0) return

  rnf=max(10.0**(-ndigit),epsm)

  analtl=max(1.0e-2,sqrt(rnf))

  if ( mod(msg/8,2)==1) go to 15
  write(ipr,901)
  write(ipr,900) (typsiz(i),i=1,n)
  write(ipr,902)
  write(ipr,900) (sx(i),i=1,n)
  write(ipr,903) fscale
  write(ipr,904) ndigit,iagflg,iahflg,iexp,method,itnlim,epsm
  write(ipr,905) stepmx,steptl,gradtl,dlt,rnf,analtl
   15 continue
!
! evaluate fcn(x)
!
  call fcn(n,x,f)
!
! evaluate analytic or finite difference gradient and check analytic
! gradient, if requested.
!
  if (iagflg /= 1) then

    call fstofd (1, 1, n, x, fcn, f, g, sx, rnf, wrk, 1)

  else

    call d1fcn (n, x, g)

    if ( mod(msg/2,2) /= 1) then
      call grdchk (n, x, fcn, f, g, typsiz, sx, fscale, rnf, analtl, wrk1, &
        msg, ipr)
      if (msg < 0) return
    end if

  end if

  call optstp(n,x,f,g,wrk1,itncnt,icscmx,itrmcd,gradtl,steptl,sx,fscale, &
    itnlim,iretcd,mxtake,ipr,msg)

  if ( itrmcd/=0) go to 700

  if ( iexp/=1) go to 80
!
! if optimization function expensive to evaluate (iexp=1), then
! hessian will be obtained by secant updates.  get initial hessian.
!
  call hsnint(nr,n,a,sx,method)
  go to 90
   80 continue
!
! evaluate analytic or finite difference hessian and check analytic
! hessian if requested (only if user-supplied analytic hessian
! routine d2fcn fills only lower triangular part and diagonal of a).
!
  if (iahflg == 1) go to 82

     if (iagflg == 1) then
       call fstofd (nr, n, n, x, d1fcn, g, a, sx, rnf, wrk1, 3)
     else
       call sndofd (nr, n, x, fcn, f, a, sx, rnf, wrk1, wrk2)
     end if

     go to 88
!
!     else
   82    if (mod(msg/4,2)==0) go to 85
!        if (mod(msg/4, 2) == 1)
!        then
        call d2fcn (nr, n, x, a)
        go to 88
!
!        else
   85     continue

          call heschk (nr, n, x, fcn, d1fcn, d2fcn, f, g, a, typsiz, &
            sx, rnf, analtl, iagflg, udiag, wrk1, wrk2, msg, ipr)
!
!           heschk evaluates d2fcn and checks it against the finite
!           difference hessian which it calculates by calling fstofd
!           (if iagflg == 1) or sndofd (otherwise).
!
        if (msg < 0) return
   88 continue

   90 continue

      if ( mod(msg/8,2)==0) then
        call result(nr,n,x,f,g,a,p,itncnt,1,ipr)
      end if
!
! iteration
!
  100 itncnt=itncnt+1
!
!  find perturbed local model hessian and its ll+ decomposition
!  (skip this step if line search or dogstep techniques being used with
!  secant updates.  cholesky decomposition l already obtained from
!  secfac.)
!
  if ( iexp==1 .and. method/=3) go to 105
  103   call chlhsn(nr,n,a,epsm,sx,udiag)
  105 continue
!
! solve for newton step:  ap=-g
!
  wrk1(1:n) = - g(1:n)

  call lltslv(nr,n,a,p,wrk1)
!
! decide whether to accept newton step  xpls=x + p
! or to choose xpls by a global strategy.
!
  if (iagflg /= 0 .or. method == 1) go to 111
  dltsav = dlt

  if ( method /= 2 ) then
    amusav = amu
    dlpsav = dltp
    phisav = phi
    phpsav = phip0
  end if

  111 continue

  if ( method==1) then
    call lnsrch(n,x,f,g,p,xpls,fpls,fcn,mxtake,iretcd, &
      stepmx,steptl,sx,ipr)
  end if

  if ( method==2) then
    call dogdrv(nr,n,x,f,g,a,p,xpls,fpls,fcn,sx,stepmx, &
      steptl,dlt,iretcd,mxtake,wrk0,wrk1,wrk2,wrk3,ipr)
  end if

  if ( method==3) then
    call hookdr(nr,n,x,f,g,a,udiag,p,xpls,fpls,fcn,sx,stepmx, &
      steptl,dlt,iretcd,mxtake,amu,dltp,phi,phip0,wrk0, &
      wrk1,wrk2,epsm,itncnt,ipr)
  end if
!
! if could not find satisfactory step and forward difference
! gradient was used, retry using central difference gradient.
!
  if ( iretcd /= 1 .or. iagflg /= 0 ) then
    go to 112
  end if
!
!  set iagflg for central differences
!
     iagflg = -1
     write(ipr,906) itncnt
     call fstocd (n, x, fcn, sx, rnf, g)
     if (method == 1) go to 105
     dlt = dltsav
     if (method == 2) go to 105
     amu = amusav
     dltp = dlpsav
     phi = phisav
     phip0 = phpsav
     go to 103
!     end if
!
! calculate step for output
!
  112 continue

  p(1:n) = xpls(1:n) - x(1:n)
!
! calculate gradient at xpls
!
  if (iagflg == (-1)) go to 116
  if (iagflg == 0) go to 118
!
! analytic gradient
  call d1fcn (n, xpls, gpls)
  go to 120
!
! central difference gradient
  116 call fstocd (n, xpls, fcn, sx, rnf, gpls)
  go to 120
!
! forward difference gradient
  118 call fstofd (1, 1, n, xpls, fcn, fpls, gpls, sx, rnf, wrk, 1)
  120 continue
!
! check whether stopping criteria satisfied
!
  call optstp(n,xpls,fpls,gpls,x,itncnt,icscmx,itrmcd,gradtl,steptl,sx,&
    fscale,itnlim,iretcd,mxtake,ipr,msg)
  if ( itrmcd/=0) go to 690
!
! evaluate hessian at xpls
!
  if ( iexp==0) go to 130

  if ( method==3) then
     call secunf(nr,n,x,g,a,udiag,xpls,gpls,epsm,itncnt,rnf, &
       iagflg,noupdt,wrk1,wrk2,wrk3)
  end if

  if ( method/=3) then
    call secfac(nr,n,x,g,a,xpls,gpls,epsm,itncnt,rnf,iagflg, &
      noupdt,wrk0,wrk1,wrk2,wrk3)
  end if

  go to 150
  130 if ( iahflg==1) go to 140

  if ( iagflg==1) then
    call fstofd(nr,n,n,xpls,d1fcn,gpls,a,sx,rnf,wrk1,3)
  else
    call sndofd(nr,n,xpls,fcn,fpls,a,sx,rnf,wrk1,wrk2)
  end if

  go to 150
  140 call d2fcn(nr,n,xpls,a)
  150 continue

  if ( mod(msg/16,2)==1) then
    call result(nr,n,xpls,fpls,gpls,a,p,itncnt,1,ipr)
  end if
!
! x <-- xpls  and  g <-- gpls  and  f <-- fpls
!
  f=fpls
  x(1:n)=xpls(1:n)
  g(1:n)=gpls(1:n)

  go to 100
!
! termination
!
! reset xpls,fpls,gpls,  if previous iterate solution
!
  690 if ( itrmcd/=3) go to 710
  700 continue

  fpls=f
  xpls(1:n)=x(1:n)
  gpls(1:n)=g(1:n)
!
! print results
!
  710 continue

  if ( mod(msg/8,2)==0) then
    call result(nr,n,xpls,fpls,gpls,a,p,itncnt,0,ipr)
  end if

  msg=0
  return
!
  900 format(' optdrv       ',5(e20.13,3x))
  901 format('0optdrv    typical x')
  902 format(' optdrv    diagonal scaling matrix for x')
  903 format(' optdrv    typical f =',e20.13)
  904 format('0optdrv    number of good digits in fcn=',i5/ &
             ' optdrv    gradient flag  =',i5,'   (=1 if analytic', &
             ' gradient supplied)'/ &
             ' optdrv    hessian flag   =',i5,'   (=1 if analytic', &
             ' hessian supplied)'/ &
             ' optdrv    expense flag   =',i5, '   (=1 if', &
             ' minimization function expensive to evaluate)'/ &
             ' optdrv    method to use  =',i5,'   (=1,2,3 for line', &
             ' search, double dogleg, more-hebdon respectively)'/ &
             ' optdrv    iteration limit=',i5/ &
             ' optdrv    machine epsilon=',e20.13)

  905 format('0optdrv    maximum step size =',e20.13/ &
             ' optdrv    step tolerance    =',e20.13/ &
             ' optdrv    gradient tolerance=',e20.13/ &
             ' optdrv    trust reg radius  =',e20.13/ &
             ' optdrv    rel noise in fcn  =',e20.13/ &
             ' optdrv    anal-fd tolerance =',e20.13)

  906 format(' optdrv    shift from forward to central differences', &
     ' in iteration ', i5)
end
subroutine optif0(nr,n,x,fcn,xpls,fpls,gpls,itrmcd,a,wrk)
!
!*******************************************************************************
!
!! OPTIF0 provides simplest interface to minimization package.
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> initial estimate of minimum
! fcn          --> name of routine to evaluate minimization function.
!                  must be declared external in calling routine.
! xpls(n)     <--  local minimum
! fpls        <--function value at local minimum xpls
! gpls(n)     <--  gradient at local minimum xpls
! itrmcd      <--  termination code
! a(n,n)       --> workspace
! wrk(n,9)     --> workspace
!
  integer n
  integer nr
!
  real a(nr,n)
  real dlt
  real fscale
  real fpls
  real gpls(n)
  real gradtl
  integer iagflg
  integer iahflg
  integer iexp
  integer ipr
  integer itnlim
  integer itrmcd
  integer method
  integer msg
  integer ndigit
  real stepmx
  real steptl
  real wrk(nr,9)
  real x(n)
  real xpls(n)
!
  external fcn,d1fcn,d2fcn
!
! equivalence wrk(n,1) = udiag(n)
!             wrk(n,2) = g(n)
!             wrk(n,3) = p(n)
!             wrk(n,4) = typsiz(n)
!             wrk(n,5) = sx(n)
!             wrk(n,6) = wrk0(n)
!             wrk(n,7) = wrk1(n)
!             wrk(n,8) = wrk2(n)
!             wrk(n,9) = wrk3(n)
!
  call dfault(n,x,wrk(1,4),fscale,method,iexp,msg,ndigit, &
    itnlim,iagflg,iahflg,ipr,dlt,gradtl,stepmx,steptl)

  call optdrv(nr,n,x,fcn,d1fcn,d2fcn,wrk(1,4),fscale, &
    method,iexp,msg,ndigit,itnlim,iagflg,iahflg,ipr, &
    dlt,gradtl,stepmx,steptl,xpls,fpls,gpls,itrmcd, &
    a,wrk(1,1),wrk(1,2),wrk(1,3),wrk(1,5),wrk(1,6), &
    wrk(1,7),wrk(1,8),wrk(1,9))

  return
end
subroutine optstp(n,xpls,fpls,gpls,x,itncnt,icscmx, &
  itrmcd,gradtl,steptl,sx,fscale,itnlim,iretcd,mxtake,ipr,msg)
!
!*******************************************************************************
!
!! OPTSTP: unconstrained minimization stopping criteria
!
!
! find whether the algorithm should terminate, due to any
! of the following:
! 1) problem solved within user tolerance
! 2) convergence within user tolerance
! 3) iteration limit reached
! 4) divergence or too restrictive maximum step (stepmx) suspected
!
!
! parameters
!
! n            --> dimension of problem
! xpls(n)      --> new iterate x[k]
! fpls         --> function value at new iterate f(xpls)
! gpls(n)      --> gradient at new iterate, g(xpls), or approximate
! x(n)         --> old iterate x[k-1]
! itncnt       --> current iteration k
! icscmx      <--> number consecutive steps >= stepmx
!                  [retain value between successive calls]
! itrmcd      <--  termination code
! gradtl       --> tolerance at which relative gradient considered close
!                  enough to zero to terminate algorithm
! steptl       --> relative step size at which successive iterates
!                  considered close enough to terminate algorithm
! sx(n)        --> diagonal scaling matrix for x
! fscale       --> estimate of scale of objective function
! itnlim       --> maximum number of allowable iterations
! iretcd       --> return code
! mxtake       --> boolean flag indicating step of maximum length used
! ipr          --> device to which to send output
! msg          --> if msg includes a term 8, suppress output
!
!
  integer n
!
  real d
  real fpls
  real fscale
  real gpls(n)
  real gradtl
  integer i
  integer icscmx
  integer ipr
  integer iretcd
  integer itncnt
  integer itnlim
  integer itrmcd
  integer jtrmcd
  integer msg
  logical mxtake
  real relgrd
  real relstp
  real rgx
  real rsx
  real steptl
  real sx(n)
  real x(n)
  real xpls(n)
!
  itrmcd=0
!
! last global step failed to locate a point lower than x.
!
  if ( iretcd == 1 ) then
    jtrmcd=3
    go to 600
  end if
!
!  Find direction in which relative gradient maximum.
!  Check whether within tolerance
!
  d=max(abs(fpls),fscale)
  rgx=0.0E+00

  do i=1,n
    relgrd=abs(gpls(i))*max(abs(xpls(i)),1./sx(i))/d
    rgx=max(rgx,relgrd)
  end do

  jtrmcd=1
  if ( rgx<=gradtl) go to 600

  if ( itncnt==0) return
!
! find direction in which relative stepsize maximum
! check whether within tolerance.
!
  rsx=0.0E+00
  do i=1,n
    relstp=abs(xpls(i)-x(i))/max(abs(xpls(i)),1./sx(i))
    rsx=max(rsx,relstp)
  end do

  jtrmcd=2
  if ( rsx<=steptl) go to 600
!
! check iteration limit
!
  jtrmcd=4
  if ( itncnt>=itnlim) go to 600
!
! check number of consecutive steps \ stepmx
!
  if ( .not.mxtake) then
    icscmx=0
    return
  else
    if (mod(msg/8,2) == 0) write(ipr,900)
    icscmx=icscmx+1
    if ( icscmx<5) return
    jtrmcd=5
  end if
!
! print termination code
!
  600 itrmcd=jtrmcd
  if (mod(msg/8,2) == 0) go to(601,602,603,604,605), itrmcd
  go to 700
  601 write(ipr,901)
  go to 700
  602 write(ipr,902)
  go to 700
  603 write(ipr,903)
  go to 700
  604 write(ipr,904)
  go to 700
  605 write(ipr,905)
!
  700 return
!
  900 format('0optstp    step of maximum length (stepmx) taken')
  901 format('0optstp    relative gradient close to zero.'/ &
             ' optstp    current iterate is probably solution.')
  902 format('0optstp    successive iterates within tolerance.'/ &
             ' optstp    current iterate is probably solution')
  903 format('0optstp    last global step failed to locate a point', &
             ' lower than x.'/ &
             ' optstp    either x is an approximate local minimum', &
             ' of the function',/ &
             ' optstp    the function is too non-linear for this algorithm,'/ &
             ' optstp    or steptl is too large.')
  904 format('optstp    iteration limit exceeded.'/'optstp    algorithm failed.')
  905 format('0optstp    maximum step size exceeded 5 consecutive times.'/ &
             ' optstp    either the function is unbounded below',/ &
             ' optstp    becomes asymptotic to a finite value', &
             ' from above in some direction',/ &
             ' optstp    or stepmx is too small')
end
subroutine passb ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSB is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  integer i
  integer idij
  integer idj
  integer idl
  integer idlj
  integer idp
  integer ik
  integer inc
  integer ipph
  integer j
  integer jc
  integer k
  integer l
  integer lc
  integer nac
  integer nt
  real wa(*)
  real wai
  real war
!
  nt = ip * idl1
  ipph = ( ip + 1 ) / 2
  idp = ip * ido

  if ( ido >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      do k = 1, l1
        ch(1:ido,k,j)  = cc(1:ido,j,k) + cc(1:ido,jc,k)
        ch(1:ido,k,jc) = cc(1:ido,j,k) - cc(1:ido,jc,k)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      do i = 1, ido
        ch(i,1:l1,j)  = cc(i,j,1:l1) + cc(i,jc,1:l1)
        ch(i,1:l1,jc) = cc(i,j,1:l1) - cc(i,jc,1:l1)
      end do
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l) = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =            wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j
      idlj = idlj + inc
      if ( idlj > idp ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  nac = 1

  if ( ido == 2 ) then
    return
  end if

  nac = 0
  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1:2,1:l1,2:ip) = ch(1:2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passb2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSB2 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,2,l1)
  real ch(ido,l1,2)
  integer i
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   + cc(i,2,k)
        ti2         = cc(i,1,k)   - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 + wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 - wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passb3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSB3 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,3,l1)
  real ch(ido,l1,3)
  real ci2
  real ci3
  real cr2
  real cr3
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer k
  real, parameter :: taui = 0.866025403784439E+00
  real, parameter :: taur = -0.5E+00
  real ti2
  real tr2
  real wa1(ido)
  real wa2(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k) - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSB4 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,4,l1)
  real ch(ido,l1,4)
  real ci1
  real ci2
  real ci3
  real ci4
  real cr1
  real cr2
  real cr3
  real cr4
  integer i
  integer k
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,4,k) - cc(2,2,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,2,k) - cc(1,4,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k) - cc(i,3,k)
        ti2 = cc(i,1,k) + cc(i,3,k)
        ti3 = cc(i,2,k) + cc(i,4,k)
        tr4 = cc(i,4,k) - cc(i,2,k)

        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,2,k) - cc(i-1,4,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3 = tr2 - tr3
        ch(i,k,1) = ti2 + ti3
        ci3 = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 - wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 + wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 - wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 + wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 - wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 + wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSB5 is a lower level routine used by CFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,5,l1)
  real ch(ido,l1,5)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer k
  real, parameter :: ti11 = 0.951056516295154E+00
  real, parameter :: ti12 = 0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real, parameter :: tr11 = 0.309016994374947E+00
  real, parameter :: tr12 = -0.809016994374947E+00
  real tr2
  real tr3
  real tr4
  real tr5
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 - wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 + wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 - wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 + wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 - wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 + wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 - wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 + wa4(i) * dr5

      end do
    end do

  end if

  return
end
subroutine passf ( nac, ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! PASSF is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  integer i
  integer idij
  integer idj
  integer idl
  integer idlj
  integer idp
  integer ik
  integer inc
  integer ipph
  integer j
  integer jc
  integer k
  integer l
  integer lc
  integer nac
  integer nt
  real wa(*)
  real wai
  real war
!
  nt = ip * idl1
  ipph = (ip+1) / 2
  idp = ip * ido

  if ( ido >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  else

    do j = 2, ipph
      jc = ip + 2 - j
      ch(1:ido,1:l1,j)  = cc(1:ido,j,1:l1) + cc(1:ido,jc,1:l1)
      ch(1:ido,1:l1,jc) = cc(1:ido,j,1:l1) - cc(1:ido,jc,1:l1)
    end do

    ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  end if

  idl = 2 - ido
  inc = 0

  do l = 2, ipph

    lc = ip + 2 - l
    idl = idl + ido

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + wa(idl-1) * ch2(ik,2)
      c2(ik,lc) =           - wa(idl)   * ch2(ik,ip)
    end do

    idlj = idl
    inc = inc + ido

    do j = 3, ipph

      jc = ip + 2 - j

      idlj = idlj + inc
      if ( idlj > idp ) then
        idlj = idlj - idp
      end if

      war = wa(idlj-1)
      wai = wa(idlj)

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + war * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) - wai * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    do ik = 2, idl1, 2
      ch2(ik-1,j)  = c2(ik-1,j) - c2(ik,jc)
      ch2(ik-1,jc) = c2(ik-1,j) + c2(ik,jc)
      ch2(ik,j)    = c2(ik,j)   + c2(ik-1,jc)
      ch2(ik,jc)   = c2(ik,j)   - c2(ik-1,jc)
    end do
  end do

  if ( ido == 2 ) then
    nac = 1
    return
  end if

  nac = 0

  c2(1:idl1,1)    = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)
  c1(2,1:l1,2:ip) = ch(2,1:l1,2:ip)

  if ( ( ido / 2 ) <= l1 ) then

    idij = 0
    do j = 2, ip
      idij = idij + 2
      do i = 4, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) + wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   - wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    idj = 2 - ido

    do j = 2, ip
      idj = idj + ido
      do k = 1, l1
        idij = idj
        do i = 4, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) + wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   - wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine passf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! PASSF2 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,2,l1)
  real ch(ido,l1,2)
  integer i
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  if ( ido <= 2 ) then

    ch(1,1:l1,1) = cc(1,1,1:l1) + cc(1,2,1:l1)
    ch(1,1:l1,2) = cc(1,1,1:l1) - cc(1,2,1:l1)
    ch(2,1:l1,1) = cc(2,1,1:l1) + cc(2,2,1:l1)
    ch(2,1:l1,2) = cc(2,1,1:l1) - cc(2,2,1:l1)

  else

    do k = 1, l1
      do i = 2, ido, 2

        ch(i-1,k,1) = cc(i-1,1,k) + cc(i-1,2,k)
        tr2         = cc(i-1,1,k) - cc(i-1,2,k)

        ch(i,k,1) = cc(i,1,k) + cc(i,2,k)
        ti2       = cc(i,1,k) - cc(i,2,k)

        ch(i,k,2)   = wa1(i-1) * ti2 - wa1(i) * tr2
        ch(i-1,k,2) = wa1(i-1) * tr2 + wa1(i) * ti2

      end do
    end do

  end if

  return
end
subroutine passf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! PASSF3 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,3,l1)
  real ch(ido,l1,3)
  real ci2
  real ci3
  real cr2
  real cr3
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer k
  real, parameter :: taui = -0.866025403784439E+00
  real, parameter :: taur = -0.5E+00
  real ti2
  real tr2
  real wa1(ido)
  real wa2(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      tr2 = cc(1,2,k) + cc(1,3,k)
      cr2 = cc(1,1,k) + taur * tr2
      ch(1,k,1) = cc(1,1,k) + tr2

      ti2 = cc(2,2,k) + cc(2,3,k)
      ci2 = cc(2,1,k) + taur * ti2
      ch(2,k,1) = cc(2,1,k) + ti2

      cr3 = taui * ( cc(1,2,k) - cc(1,3,k) )
      ci3 = taui * ( cc(2,2,k) - cc(2,3,k) )

      ch(1,k,2) = cr2 - ci3
      ch(1,k,3) = cr2 + ci3
      ch(2,k,2) = ci2 + cr3
      ch(2,k,3) = ci2 - cr3

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        tr2 = cc(i-1,2,k) + cc(i-1,3,k)
        cr2 = cc(i-1,1,k) + taur * tr2
        ch(i-1,k,1) = cc(i-1,1,k) + tr2

        ti2 = cc(i,2,k) + cc(i,3,k)
        ci2 = cc(i,1,k) + taur * ti2
        ch(i,k,1) = cc(i,1,k) + ti2

        cr3 = taui * ( cc(i-1,2,k) - cc(i-1,3,k) )
        ci3 = taui * ( cc(i,2,k)   - cc(i,3,k) )

        dr2 = cr2 - ci3
        dr3 = cr2 + ci3
        di2 = ci2 + cr3
        di3 = ci2 - cr3

        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3

      end do
    end do

  end if

  return
end
subroutine passf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! PASSF4 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,4,l1)
  real ch(ido,l1,4)
  real ci1
  real ci2
  real ci3
  real ci4
  real cr1
  real cr2
  real cr3
  real cr4
  integer i
  integer k
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti1 = cc(2,1,k) - cc(2,3,k)
      ti2 = cc(2,1,k) + cc(2,3,k)
      tr4 = cc(2,2,k) - cc(2,4,k)
      ti3 = cc(2,2,k) + cc(2,4,k)
      tr1 = cc(1,1,k) - cc(1,3,k)
      tr2 = cc(1,1,k) + cc(1,3,k)
      ti4 = cc(1,4,k) - cc(1,2,k)
      tr3 = cc(1,2,k) + cc(1,4,k)

      ch(1,k,1) = tr2 + tr3
      ch(1,k,3) = tr2 - tr3
      ch(2,k,1) = ti2 + ti3
      ch(2,k,3) = ti2 - ti3
      ch(1,k,2) = tr1 + tr4
      ch(1,k,4) = tr1 - tr4
      ch(2,k,2) = ti1 + ti4
      ch(2,k,4) = ti1 - ti4

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti1 = cc(i,1,k)   - cc(i,3,k)
        ti2 = cc(i,1,k)   + cc(i,3,k)
        ti3 = cc(i,2,k)   + cc(i,4,k)
        tr4 = cc(i,2,k)   - cc(i,4,k)
        tr1 = cc(i-1,1,k) - cc(i-1,3,k)
        tr2 = cc(i-1,1,k) + cc(i-1,3,k)
        ti4 = cc(i-1,4,k) - cc(i-1,2,k)
        tr3 = cc(i-1,2,k) + cc(i-1,4,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 + tr4
        cr4 = tr1 - tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-1) * cr2 + wa1(i) * ci2
        ch(i,k,2)   = wa1(i-1) * ci2 - wa1(i) * cr2
        ch(i-1,k,3) = wa2(i-1) * cr3 + wa2(i) * ci3
        ch(i,k,3)   = wa2(i-1) * ci3 - wa2(i) * cr3
        ch(i-1,k,4) = wa3(i-1) * cr4 + wa3(i) * ci4
        ch(i,k,4)   = wa3(i-1) * ci4 - wa3(i) * cr4

      end do
    end do

  end if

  return
end
subroutine passf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! PASSF5 is a lower level routine used by CFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,5,l1)
  real ch(ido,l1,5)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer k
  real, parameter :: ti11 = -0.951056516295154E+00
  real, parameter :: ti12 = -0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real tr2
  real tr3
  real tr4
  real tr5
  real, parameter :: tr11 =  0.309016994374947E+00
  real, parameter :: tr12 = -0.809016994374947E+00
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  if ( ido == 2 ) then

    do k = 1, l1

      ti5 = cc(2,2,k) - cc(2,5,k)
      ti2 = cc(2,2,k) + cc(2,5,k)
      ti4 = cc(2,3,k) - cc(2,4,k)
      ti3 = cc(2,3,k) + cc(2,4,k)
      tr5 = cc(1,2,k) - cc(1,5,k)
      tr2 = cc(1,2,k) + cc(1,5,k)
      tr4 = cc(1,3,k) - cc(1,4,k)
      tr3 = cc(1,3,k) + cc(1,4,k)

      ch(1,k,1) = cc(1,1,k) + tr2 + tr3
      ch(2,k,1) = cc(2,1,k) + ti2 + ti3

      cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(2,1,k) + tr11 * ti2 + tr12 * ti3
      cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(2,1,k) + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      ch(1,k,2) = cr2 - ci5
      ch(1,k,5) = cr2 + ci5
      ch(2,k,2) = ci2 + cr5
      ch(2,k,3) = ci3 + cr4
      ch(1,k,3) = cr3 - ci4
      ch(1,k,4) = cr3 + ci4
      ch(2,k,4) = ci3 - cr4
      ch(2,k,5) = ci2 - cr5

    end do

  else

    do k = 1, l1
      do i = 2, ido, 2

        ti5 = cc(i,2,k) - cc(i,5,k)
        ti2 = cc(i,2,k) + cc(i,5,k)
        ti4 = cc(i,3,k) - cc(i,4,k)
        ti3 = cc(i,3,k) + cc(i,4,k)

        tr5 = cc(i-1,2,k) - cc(i-1,5,k)
        tr2 = cc(i-1,2,k) + cc(i-1,5,k)
        tr4 = cc(i-1,3,k) - cc(i-1,4,k)
        tr3 = cc(i-1,3,k) + cc(i-1,4,k)

        ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
        ch(i,k,1)   = cc(i,1,k)   + ti2 + ti3

        cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
        ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
        cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
        ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

        cr5 = ti11 * tr5 + ti12 * tr4
        ci5 = ti11 * ti5 + ti12 * ti4
        cr4 = ti12 * tr5 - ti11 * tr4
        ci4 = ti12 * ti5 - ti11 * ti4

        dr3 = cr3 - ci4
        dr4 = cr3 + ci4
        di3 = ci3 + cr4
        di4 = ci3 - cr4
        dr5 = cr2 + ci5
        dr2 = cr2 - ci5
        di5 = ci2 - cr5
        di2 = ci2 + cr5

        ch(i-1,k,2) = wa1(i-1) * dr2 + wa1(i) * di2
        ch(i,k,2)   = wa1(i-1) * di2 - wa1(i) * dr2
        ch(i-1,k,3) = wa2(i-1) * dr3 + wa2(i) * di3
        ch(i,k,3)   = wa2(i-1) * di3 - wa2(i) * dr3
        ch(i-1,k,4) = wa3(i-1) * dr4 + wa3(i) * di4
        ch(i,k,4)   = wa3(i-1) * di4 - wa3(i) * dr4
        ch(i-1,k,5) = wa4(i-1) * dr5 + wa4(i) * di5
        ch(i,k,5)   = wa4(i-1) * di5 - wa4(i) * dr5

      end do
    end do

  end if

  return
end
subroutine pchce(ic,vc,n,x,h,slope,d,incfd,ierr)
!
!*******************************************************************************
!
!! PCHCE is called by pchic to set end derivatives as requested by the user.
!
!
!    it must be called after interior derivative values have been set.
!
!
!    to facilitate two-dimensional applications, includes an increment
!    between successive values of the d-array.
!
!  Parameters:
!
!     ic -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           ic(1) = ibeg, desired condition at beginning of data.
!           ic(2) = iend, desired condition at end of data.
!           ( see prologue to pchic for details. )
!
!     vc -- (input) real array of length 2 specifying desired boundary
!           values.  vc(1) need be set only if ic(1) = 2 or 3 .
!                    vc(2) need be set only if ic(2) = 2 or 3 .
!
!     n -- (input) number of data points.  (assumes n>=2)
!
!     x -- (input) real array of independent variable values.  (the
!           elements of x are assumed to be strictly increasing.)
!
!     h -- (input) real array of interval lengths.
!     slope -- (input) real array of data slopes.
!           if the data are (x(i),y(i)), i=1(1)n, then these inputs are:
!                  h(i) =  x(i+1)-x(i),
!              slope(i) = (y(i+1)-y(i))/h(i),  i=1(1)n-1.
!
!     d -- (input) real array of derivative values at the data points.
!           the value corresponding to x(i) must be stored in
!                d(1+(i-1)*incfd),  i=1(1)n.
!          (output) the value of d at x(1) and/or x(n) is changed, if
!           necessary, to produce the requested boundary conditions.
!           no other entries in d are changed.
!
!     incfd -- (input) increment between successive values in d.
!           this argument is provided primarily for 2-d applications.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning errors:
!              ierr = 1  if ibeg<0 and d(1) had to be adjusted for
!                        monotonicity.
!              ierr = 2  if iend<0 and d(1+(n-1)*incfd) had to be
!                        adjusted for monotonicity.
!              ierr = 3  if both of the above are true.
!
!  programming notes:
!
!     1. the function  pchst(arg1,arg2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!     3. one could reduce the number of arguments and amount of local
!        storage, at the expense of reduced code clarity, by passing in
!        the array wk (rather than splitting it into h and slope) and
!        increasing its length enough to incorporate stemp and xtemp.
!     4. the two monotonicity checks only use the sufficient conditions.
!        thus, it is possible (but unlikely) for a boundary condition to
!        be changed, even though the original interpolant was monotonic.
!        (at least the result is a continuous function of the data.)
!
  integer incfd
  integer n
!
  real d(incfd,n)
  real h(n)
  integer ibeg
  integer ic(2)
  integer iend
  integer ierf
  integer ierr
  integer index
  integer j
  integer k
  real pchdf
  real pchst
  real slope(n)
  real stemp(3)
  real vc(2)
  real x(n)
  real xtemp(4)
!
  ibeg = ic(1)
  iend = ic(2)
  ierr = 0
!
!  set to default boundary conditions if n is too small.
!
  if ( abs(ibeg)>n )  ibeg = 0
  if ( abs(iend)>n )  iend = 0
!
!  treat beginning boundary condition.
!
  if (ibeg == 0)  go to 2000

  k = abs(ibeg)
  if (k == 1)  then
!        boundary value provided.
     d(1,1) = vc(1)
  else if (k == 2)  then
!        boundary second derivative provided.
!
     d(1,1) = 0.5 *( (3.0*slope(1) - d(1,2)) - 0.5*vc(1)*h(1) )
  else if (k < 5)  then
!        use k-point derivative formula.
!        pick up first k points, in reverse order.
!
     do j = 1, k
        index = k-j+1
        xtemp(j) = x(index)
        if (j < k)  stemp(j) = slope(index-1)
     end do

     d(1,1) = pchdf (k, xtemp, stemp, ierf)

     if (ierf /= 0)  go to 5001
  else
!        use 'not a knot' condition.
     d(1,1) = ( 3.0*(h(1)*slope(2) + h(2)*slope(1)) &
       - 2.0E+00 * (h(1)+h(2))*d(1,2) - h(1)*d(1,3) ) / h(2)
  end if
!
!  check d(1,1) for compatibility with monotonicity.
!
  if ( ibeg <= 0 ) then

    if (slope(1) == 0.0)  then
      if (d(1,1) /= 0.0)  then
        d(1,1) = 0.0E+00
        ierr = ierr + 1
      end if
    else if ( pchst(d(1,1),slope(1)) < 0.0)  then
      d(1,1) = 0.0E+00
      ierr = ierr + 1
    else if ( abs(d(1,1)) > 3.0*abs(slope(1)) )  then
      d(1,1) = 3.0*slope(1)
      ierr = ierr + 1
    end if

  end if

2000 continue
!
!  treat end boundary condition.
!
  if (iend == 0)  go to 5000
  k = abs(iend)
  if (k == 1)  then
!        boundary value provided.
     d(1,n) = vc(2)
  else if (k == 2)  then
!        boundary second derivative provided.
     d(1,n) = 0.5 *( (3.0*slope(n-1) - d(1,n-1)) + 0.5 *vc(2)*h(n-1) )
  else if (k < 5)  then
!
!  use k-point derivative formula.  pick up last k points.
!
     do j = 1, k
        index = n-k+j
        xtemp(j) = x(index)
        if (j < k)  stemp(j) = slope(index)
     end do

     d(1,n) = pchdf (k, xtemp, stemp, ierf)

     if (ierf /= 0)  go to 5001
  else
!        use 'not a knot' condition.
     d(1,n) = ( 3.0*(h(n-1)*slope(n-2) + h(n-2)*slope(n-1)) &
       - 2.0E+00 * (h(n-1)+h(n-2))*d(1,n-1) - h(n-1)*d(1,n-2) ) / h(n-2)
  end if
!
  if (iend > 0)  go to 5000
!
!  check d(1,n) for compatibility with monotonicity.
!
  if (slope(n-1) == 0.0)  then
     if (d(1,n) /= 0.0)  then
        d(1,n) = 0.0E+00
        ierr = ierr + 2
     end if
  else if ( pchst(d(1,n),slope(n-1)) < 0.0)  then
     d(1,n) = 0.0E+00
     ierr = ierr + 2
  else if ( abs(d(1,n)) > 3.0*abs(slope(n-1)) )  then
     d(1,n) = 3.0*slope(n-1)
     ierr = ierr + 2
  end if
!
!  normal return.
!
 5000 continue
  return
!
!  error return.
!
 5001 continue
!     error return from pchdf.
!  This case should never occur.
!
  ierr = -1
  call xerror ('pchce -- error return from pchdf', 32, ierr, 1)
  return
end
subroutine pchci(n,h,slope,d,incfd)
!
!*******************************************************************************
!
!! PCHCI sets derivatives for a monotone piecewise cubic hermite interpolant.
!
!
!  Discussion:
!
!    default boundary conditions are provided which are compatible
!    with monotonicity.  if the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    to facilitate two-dimensional applications, includes an increment
!    between successive values of the d-array.
!
!    the resulting piecewise cubic hermite function should be identical
!    (within roundoff error) to that produced by pchim.
!
!  Parameters:
!
!     n -- (input) number of data points.
!           if n=2, simply does linear interpolation.
!
!     h -- (input) real array of interval lengths.
!     slope -- (input) real array of data slopes.
!           if the data are (x(i),y(i)), i=1(1)n, then these inputs are:
!                  h(i) =  x(i+1)-x(i),
!              slope(i) = (y(i+1)-y(i))/h(i),  i=1(1)n-1.
!
!     d -- (output) real array of derivative values at the data points.
!           if the data are monotonic, these values will determine a
!           a monotone cubic hermite function.
!           the value corresponding to x(i) is stored in
!                d(1+(i-1)*incfd),  i=1(1)n.
!           no other entries in d are changed.
!
!     incfd -- (input) increment between successive values in d.
!           this argument is provided primarily for 2-d applications.
!
!  programmed by:  fred n. fritsch,  fts 532-4275, (415) 422-4275,
!                  mathematics and statistics division,
!                  lawrence livermore national laboratory.
!
!  programming notes:
!
!     1. the function  pchst(arg1,arg2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!
  integer incfd
  integer n
!
  real d(incfd,n)
  real del1
  real del2
  real dmax
  real dmin
  real drat1
  real drat2
  real h(n)
  real hsum
  real hsumt3
  integer i
  integer nless1
  real pchst
  real slope(n)
  real w1
  real w2
!
  nless1 = n - 1
  del1 = slope(1)
!
!  special case n=2 -- use linear interpolation.
!
  if (nless1 > 1)  go to 10
  d(1,1) = del1
  d(1,n) = del1
  go to 5000
!
!  normal case  (n >= 3).
!
   10 continue
  del2 = slope(2)
!
!  set d(1) via non-centered three-point formula, adjusted to be shape-preserving.
!
  hsum = h(1) + h(2)
  w1 = (h(1) + hsum)/hsum
  w2 = -h(1)/hsum
  d(1,1) = w1*del1 + w2*del2

  if ( pchst(d(1,1),del1) <= 0.0)  then
     d(1,1) = 0.0E+00
  else if ( pchst(del1,del2) < 0.0)  then
!        need do this check only if monotonicity switches.
     dmax = 3.0*del1
     if (abs(d(1,1)) > abs(dmax))  d(1,1) = dmax
  end if
!
!  loop through interior points.
!
  do i = 2, nless1

    if (i /= 2)  then
      hsum = h(i-1) + h(i)
      del1 = del2
      del2 = slope(i)
    end if
!
!  set d(i)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0E+00
!
!  use brodlie modification of butland formula.
!
    if ( pchst(del1,del2) > 0.0E+00 ) then

      hsumt3 = hsum+hsum+hsum
      w1 = (hsum + h(i-1))/hsumt3
      w2 = (hsum + h(i)  )/hsumt3
      dmax = max( abs(del1), abs(del2) )
      dmin = min ( abs(del1), abs(del2) )
      drat1 = del1/dmax
      drat2 = del2/dmax
      d(1,i) = dmin/(w1*drat1 + w2*drat2)

    end if

  end do
!
!  set d(n) via non-centered three-point formula, adjusted to be shape-preserving.
!
  w1 = -h(n-1)/hsum
  w2 = (h(n-1) + hsum)/hsum
  d(1,n) = w1*del1 + w2*del2

  if ( pchst(d(1,n),del2) <= 0.0)  then
    d(1,n) = 0.0E+00
  else if ( pchst(del1,del2) < 0.0)  then
    dmax = 3.0*del2
    if (abs(d(1,n)) > abs(dmax)) then
      d(1,n) = dmax
    end if
  end if

 5000 continue
  return
end
subroutine pchcs(switch,n,h,slope,d,incfd,ierr)
!
!*******************************************************************************
!
!! PCHCS is called by  pchic  to adjust the values of d in the vicinity of a
!     switch in direction of monotonicity, to produce a more "visually
!     pleasing" curve than that given by  pchim .
!
!  Parameters:
!
!     switch -- (input) indicates the amount of control desired over
!           local excursions from data.
!
!     n -- (input) number of data points.  (assumes n>2 .)
!
!     h -- (input) real array of interval lengths.
!     slope -- (input) real array of data slopes.
!           if the data are (x(i),y(i)), i=1(1)n, then these inputs are:
!                  h(i) =  x(i+1)-x(i),
!              slope(i) = (y(i+1)-y(i))/h(i),  i=1(1)n-1.
!
!     d -- (input) real array of derivative values at the data points,
!           as determined by pchci.
!          (output) derivatives in the vicinity of switches in direction
!           of monotonicity may be adjusted to produce a more "visually
!           pleasing" curve.
!           the value corresponding to x(i) is stored in
!                d(1+(i-1)*incfd),  i=1(1)n.
!           no other entries in d are changed.
!
!     incfd -- (input) increment between successive values in d.
!           this argument is provided primarily for 2-d applications.
!
!     ierr -- (output) error flag.  should be zero.
!           if negative, trouble in pchsw.  (should never happen.)
!
!    warning:  this routine does no validity-checking of arguments.
!
!  programmed by:  fred n. fritsch,  fts 532-4275, (415) 422-4275,
!                  mathematics and statistics division,
!                  lawrence livermore national laboratory.
!
!  programming notes:
!
!     1. the function  pchst(arg1,arg2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!
  integer incfd
  integer n
!
  real d(incfd,n)
  real h(n)
  integer ierr
  real slope(n)
  real switch

  integer  i, indx, k, nless1
  real  del(3), dext, dfloc, dfmx, fact, fudge, slmax
  real wtave(2)
  real  pchst
!
!  define inline function for weighted average of slopes.
!
  real  pchsd, s1, s2, h1, h2
  pchsd(s1,s2,h1,h2) = (h2/(h1+h2))*s1 + (h1/(h1+h2))*s2
!
!  initialize.
!
  data  fudge /4./
!
  ierr = 0
  nless1 = n - 1
!
!  loop over segments.
!
  do 900  i = 2, nless1

     if ( pchst(slope(i-1),slope(i)) )  100, 300, 900

  100    continue
!
!  slope switches monotonicity at i-th point
!
!           do not change d if 'up-down-up'.
        if (i > 2)  then
           if ( pchst(slope(i-2),slope(i)) > 0.0)  go to 900

        end if
        if (i < nless1)  then
           if ( pchst(slope(i+1),slope(i-1)) > 0.0)  go to 900

        end if
!
!  compute provisional value for d(1,i).
!
        dext = pchsd (slope(i-1), slope(i), h(i-1), h(i))
!
!  determine which interval contains the extremum.
!
        if ( pchst(dext, slope(i-1)) )  200, 900, 250
!
  200       continue
!              dext and slope(i-1) have opposite signs --
!                        extremum is in (x(i-1),x(i)).
           k = i-1
!              set up to compute new values for d(1,i-1) and d(1,i).
           wtave(2) = dext
           if (k > 1) then
             wtave(1) = pchsd (slope(k-1), slope(k), h(k-1), h(k))
           end if
           go to 400
!
  250       continue
!              dext and slope(i) have opposite signs --
!                        extremum is in (x(i),x(i+1)).
           k = i
!              set up to compute new values for d(1,i) and d(1,i+1).
           wtave(1) = dext
           if (k < nless1) then
             wtave(2) = pchsd (slope(k), slope(k+1), h(k), h(k+1))
           end if
           go to 400
!
  300    continue
!
!  at least one of slope(i-1) and slope(i) is zero
!  check for flat-topped peak
!
        if (i == nless1)  go to 900
        if ( pchst(slope(i-1), slope(i+1)) >= 0.0)  go to 900

!
!           we have flat-topped peak on (x(i),x(i+1)).
!
        k = i
!           set up to compute new values for d(1,i) and d(1,i+1).
        wtave(1) = pchsd (slope(k-1), slope(k), h(k-1), h(k))
        wtave(2) = pchsd (slope(k), slope(k+1), h(k), h(k+1))
!
  400    continue
!
!  at this point we have determined that there will be an extremum
!        on (x(k),x(k+1)), where k=i or i-1, and have set array wtave--
!           wtave(1) is a weighted average of slope(k-1) and slope(k),
!                    if k>1
!           wtave(2) is a weighted average of slope(k) and slope(k+1),
!                    if k<n-1
!
     slmax = abs(slope(k))
     if (k > 1)    slmax = max ( slmax, abs(slope(k-1)) )
     if (k<nless1) slmax = max ( slmax, abs(slope(k+1)) )
!
     if (k > 1)  del(1) = slope(k-1) / slmax
     del(2) = slope(k) / slmax
     if (k<nless1)  del(3) = slope(k+1) / slmax
!
     if ((k>1) .and. (k<nless1))  then
!           normal case -- extremum is not in a boundary interval.
        fact = fudge* abs(del(3)*(del(1)-del(2))*(wtave(2)/slmax))
        d(1,k) = d(1,k) + min (fact,1.0)*(wtave(1) - d(1,k))
        fact = fudge* abs(del(1)*(del(3)-del(2))*(wtave(1)/slmax))
        d(1,k+1) = d(1,k+1) + min (fact,1.0)*(wtave(2) - d(1,k+1))
     else
!           special case k=1 (which can occur only if i=2) or
!                        k=nless1 (which can occur only if i=nless1).
        fact = fudge* abs(del(2))
        d(1,i) = min (fact,1.0) * wtave(i-k+1)
!              note that i-k+1 = 1 if k=i  (=nless1),
!                        i-k+1 = 2 if k=i-1(=1).
     end if
!
!  adjust if necessary to limit excursions from data.
!
     if (switch <= 0.0)  go to 900

     dfloc = h(k)*abs(slope(k))
     if (k > 1)    dfloc = max ( dfloc, h(k-1)*abs(slope(k-1)) )
     if (k<nless1) dfloc = max ( dfloc, h(k+1)*abs(slope(k+1)) )
     dfmx = switch*dfloc
     indx = i-k+1
!        indx = 1 if k=i, 2 if k=i-1.

     call pchsw (dfmx, indx, d(1,k), d(1,k+1), h(k), slope(k), ierr)

     if (ierr /= 0)  return

  900 continue

  return
end
function pchdf ( k, x, s, ierr )
!
!*******************************************************************************
!
!! PCHDF approximates a derivative using divided differences of data.
!
!
!  Discussion:
!
!    The routine uses a divided difference formulation to compute a K-point
!    approximation to the derivative at X(K) based on the data in X and S.
!
!    It is called by PCHCE and PCHSP to compute 3- and 4-point boundary
!    derivative approximations.
!
!  Reference:
!
!    Carl de Boor,
!    a practical guide to splines,
!    springer-verlag (new york, 1978), pp. 10-16.
!
!  Parameters:
!
!     on input:
!        k      is the order of the desired derivative approximation.
!               k must be at least 3 (error return if not).
!        x      contains the k values of the independent variable.
!               x need not be ordered, but the values **must** be
!               distinct.  (not checked here.)
!        s      contains the associated slope values:
!                  s(i) = (f(i+1)-f(i))/(x(i+1)-x(i)), i=1(1)k-1.
!               (note that s need only be of length k-1.)
!
!     on return:
!        s      will be destroyed.
!        ierr   will be set to -1 if k<2 .
!        pchdf  will be set to the desired derivative approximation if
!               ierr=0 or to zero if ierr=-1.
!
!  programmed by:  fred n. fritsch,  fts 532-4275, (415) 422-4275,
!                  mathematics and statistics division,
!                  lawrence livermore national laboratory.
!
  integer k
!
  integer i
  integer ierr
  integer j
  real pchdf
  real s(k)
  real value
  real x(k)
!
!  check for legal value of k.
!
  if (k < 3)  go to 5001
!
!  compute coefficients of interpolating polynomial.
!
  do j = 2, k-1
    do i = 1, k-j
      s(i) = (s(i+1)-s(i))/(x(i+j)-x(i))
    end do
  end do
!
!  evaluate derivative at x(k).
!
  value = s(1)

  do i = 2, k-1
    value = s(i) + value * ( x(k) - x(i) )
  end do
!
!  normal return.
!
  ierr = 0
  pchdf = value
  return
!
!  error return.
!
 5001 continue
!     k<3 return.
  ierr = -1
  call xerror ('pchdf -- k less than three', 26, ierr, 1)
  pchdf = 0.0E+00
  return
end
subroutine pchev ( n, x, f, d, nval, xval, fval, dval, ierr )
!
!*******************************************************************************
!
!! PCHEV evaluates the function and first derivative of a piecewise
!            cubic hermite or spline function at an array of points xval,
!            easy to use.
!
!  Description
!
!          pchev:  piecewise cubic hermite or spline derivative evaluator,
!                  easy to use.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!     evaluates the function and first derivative of the cubic hermite
!     or spline function defined by  n, x, f, d, at the array of points xval.
!
!     this is an easy to use driver for the routines by f.n. fritsch
!     described in reference (2) below. those also have other capabilities.
!
!  Parameters:
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!             x(i-1) < x(i),  i = 2(1)n. (error return if not.)
!
!     f -- (input) real array of function values.  f(i) is
!           the value corresponding to x(i).
!
!     d -- (input) real array of derivative values.  d(i) is
!           the value corresponding to x(i).
!
!  nval -- (input) number of points at which the functions are to be
!           evaluated. ( error return if nval<1 )
!
!  xval -- (input) real array of points at which the functions are to
!           be evaluated.
!
!          notes:
!           1. the evaluation will be most efficient if the elements
!              of xval are increasing relative to x;
!              that is,   xval(j) >= x(i)
!              implies    xval(k) >= x(i),  all k>=j .
!           2. if any of the xval are outside the interval [x(1),x(n)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!  fval -- (output) real array of values of the cubic hermite function
!           defined by  n, x, f, d  at the points  xval.
!
!  dval -- (output) real array of values of the first derivative of
!           the same function at the points  xval.
!
!  ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning error:
!              ierr>0  means that extrapolation was performed at
!                 ierr points.
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -3  if the x-array is not strictly increasing.
!              ierr = -4  if nval<1 .
!           (output arrays have not been changed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!              ierr = -5  if an error has occurred in the lower-level
!                         routine chfdv.  nb: this should never happen.
!                         notify the author **immediately** if it does.
!
!
!  references  1. f.n.fritsch and r.e.carlson, 'monotone piecewise
!                 cubic interpolation,' siam j.numer.anal. 17, 2 (april
!                 1980), 238-246.
!               2. f.n.fritsch, 'piecewise cubic hermite interpolation
!                 package, final specifications', lawrence livermore
!                 national laboratory, computer documentation ucid-30194,
!                 august 1982.
!
  integer n
  integer nval
!
  real d(n)
  real dval(nval)
  real f(n)
  real fval(nval)
  integer ierr
  integer, save :: incfd = 1
  logical, save :: skip = .true.
  real x(n)
  real xval(nval)
!
  call pchfd ( n, x, f, d, incfd, skip, nval, xval, fval, dval, ierr )

  return
end
subroutine pchez ( n, x, f, d, spline, wk, lwk, ierr )
!
!*******************************************************************************
!
!! PCHEZ carries out easy to use spline or cubic hermite interpolation.
!
!
!  Description
!
!          pchez:  piecewise cubic interpolation, easy to use.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!     sets derivatives for spline (two continuous derivatives) or
!     hermite cubic (one continuous derivative) interpolation.
!     spline interpolation is smoother, but may not "look" right if the
!     data contains both "steep" and "flat" sections.  hermite cubics
!     can produce a "visually pleasing" and monotone interpolant to
!     monotone data. this is an easy to use driver for the routines
!     by f. n. fritsch in reference (4) below. various boundary
!     conditions are set to default values by pchez. many other choices
!     are available in the subroutines pchic, pchim and pchsp.
!
!     use pchev to evaluate the resulting function and its derivative.
!
!  Parameters:
!
!     n -- (input) number of data points.  (error return if n<2 .)
!           if n=2, simply does linear interpolation.
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of dependent variable values to be inter-
!           polated.  f(i) is value corresponding to x(i).
!
!     d -- (output) real array of derivative values at the data points.
!
!     spline -- (input) logical variable to specify if the interpolant
!           is to be a spline with two continuous derivaties
!           (set spline=.true.) or a hermite cubic interpolant with one
!           continuous derivative (set spline=.false.).
!        note: if spline=.true. the interpolating spline satisfies the
!           default "not-a-knot" boundary condition, with a continuous
!           third derivative at x(2) and x(n-1). see reference (3).
!              if spline=.false. the interpolating hermite cubic will be
!           monotone if the input data is monotone. boundary conditions are
!           computed from the derivative of a local quadratic unless this
!           alters monotonicity.
!
!     wk -- (scratch) real work array, which must be declared by the calling
!           program to be at least 2*n if spline is .true. and not used
!           otherwise.
!
!     lwk -- (input) length of work array wk. (error return if
!           lwk<2*n and spline is .true., not checked otherwise.)
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning error:
!              ierr>0  (can only occur when spline=.false.) means that
!                 ierr switches in the direction of monotonicity were detected.
!                 when spline=.false.,  pchez guarantees that if the input
!                 data is monotone, the interpolant will be too. this warning
!                 is to alert you to the fact that the input data was not
!                 monotone.
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -3  if the x-array is not strictly increasing.
!              ierr = -7  if lwk is less than 2*n and spline is .true.
!             (the d-array has not been changed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
!  references  1. f.n.fritsch and r.e.carlson, 'monotone piecewise
!                 cubic interpolation,' siam j.numer.anal. 17, 2 (april
!                 1980), 238-246.
!               2. f.n.fritsch and j.butland, 'a method for constructing
!                 local monotone piecewise cubic interpolants,' llnl
!                 preprint ucrl-87559 (april 1982).
!               3. carl de boor, a practical guide to splines, springer-
!                 verlag (new york, 1978).  (esp. chapter iv, pp.49-62.)
!               4. f.n.fritsch, 'piecewise cubic hermite interpolation
!                 package, final specifications', lawrence livermore
!                 national laboratory, computer documentation ucid-30194,
!                 august 1982.
!
  integer lwk
  integer n
!
  real d(n)
  real f(n)
  integer, save, dimension ( 2 ) :: ic = (/ 0, 0 /)
  integer ierr
  integer, parameter :: incfd = 1
  logical spline
  real vc(2)
  real wk(lwk)
  real x(n)
!
  if ( spline ) then
    call  pchsp ( ic, vc, n, x, f, d, incfd, wk, lwk, ierr )
  else
    call  pchim ( n, x, f, d, incfd, ierr )
  end if

  return
end
subroutine pchfd(n,x,f,d,incfd,skip,ne,xe,fe,de,ierr)
!
!*******************************************************************************
!
!! PCHFD evaluates a piecewise cubic hermite function and its first
!            derivative at an array of points.  may be used by itself
!            for hermite interpolation, or as an evaluator for pchim
!            or pchic.  if only function values are required, use
!            PCHFE instead.
!
!  Description:
!
!          pchfd:  piecewise cubic hermite function and derivative
!                  evaluator
!
!     evaluates the cubic hermite function defined by  n, x, f, d,  to-
!     gether with its first derivative, at the points  xe(j), j=1(1)ne.
!
!     if only function values are required, use PCHFE, instead.
!
!     to provide compatibility with pchim and pchic, includes an
!     increment between successive values of the f- and d-arrays.
!
!  Parameters:
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of function values.  f(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     d -- (input) real array of derivative values.  d(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     incfd -- (input) increment between successive values in f and d.
!           (error return if  incfd<1 .)
!
!     skip -- (input/output) logical variable which should be set to
!           .true. if the user wishes to skip checks for validity of
!           preceding parameters, or to .false. otherwise.
!           this will save time in case these checks have already
!           been performed (say, in pchim or pchic).
!           skip will be set to .true. on normal return.
!
!     ne -- (input) number of evaluation points.  (error return if
!           ne<1 .)
!
!     xe -- (input) real array of points at which the functions are to
!           be evaluated.
!
!
!          notes:
!           1. the evaluation will be most efficient if the elements
!              of xe are increasing relative to x;
!              that is,   xe(j) >= x(i)
!              implies    xe(k) >= x(i),  all k>=j .
!           2. if any of the xe are outside the interval [x(1),x(n)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     fe -- (output) real array of values of the cubic hermite function
!           defined by  n, x, f, d  at the points  xe.
!
!     de -- (output) real array of values of the first derivative of
!           the same function at the points  xe.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning error:
!              ierr>0  means that extrapolation was performed at
!                 ierr points.
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!              ierr = -4  if ne<1 .
!           (output arrays have not been changed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!              ierr = -5  if an error has occurred in the lower-level
!                         routine chfdv.  nb: this should never happen.
!                         notify the author **immediately** if it does.
!
!  programming notes:
!
!     2. most of the coding between the call to chfdv and the end of
!        the ir-loop could be eliminated if it were permissible to
!        assume that xe is ordered relative to x.
!
!     3. chfdv does not assume that x1 is less than x2.  thus, it would
!        be possible to write a version of pchfd that assumes a strict-
!        ly decreasing x-array by simply running the ir-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. the present code has a minor bug, which i have decided is not
!        worth the effort that would be required to fix it.
!        if xe contains points in [x(n-1),x(n)], followed by points <
!        x(n-1), followed by points >x(n), the extrapolation points
!        will be counted (at least) twice in the total returned in ierr.
!
  integer incfd
  integer n
  integer ne
!
  real d(incfd,n)
  real de(ne)
  real f(incfd,n)
  real fe(ne)
  integer i
  integer ierc
  integer ierr
  integer ir
  integer j
  integer jfirst
  integer next(2)
  integer nj
  logical skip
  real x(n)
  real xe(ne)
!
!  validity-check arguments.
!
  if (skip)  go to 5

  if ( n<2 )  go to 5001
  if ( incfd<1 )  go to 5002

  do i = 2, n
    if ( x(i)<=x(i-1) ) then
      go to 5003
    end if
  end do
!
!  function definition is ok, go on.
!
    5 continue
  if ( ne<1 )  go to 5004
  ierr = 0
  skip = .true.
!
!  loop over intervals.        (   interval index is  il = ir-1  . )
!                              ( interval is x(il)<=x<x(ir) . )
  jfirst = 1
  ir = 2
   10 continue
!
!     skip out of loop if have processed all evaluation points.
!
     if (jfirst > ne)  go to 5000
!
!     locate all points in interval.
!
     do j = jfirst, ne
       if (xe(j) >= x(ir))  go to 30
     end do

     j = ne + 1
     go to 40
!
!     have located first point beyond interval.
!
   30    continue
     if (ir == n)  j = ne + 1
!
   40    continue
     nj = j - jfirst
!
!     skip evaluation if no points in interval.
!
     if (nj == 0)  go to 50
!
!     evaluate cubic at xe(i),  i = jfirst (1) j-1 .
!
    call chfdv (x(ir-1),x(ir), f(1,ir-1),f(1,ir), d(1,ir-1),d(1,ir), &
      nj, xe(jfirst), fe(jfirst), de(jfirst), next, ierc)

     if (ierc < 0)  go to 5005

     if (next(2) == 0)  go to 42
!        if (next(2) > 0)  then
!           in the current set of xe-points, there are next(2) to the
!           right of x(ir).
!
        if (ir < n)  go to 41
!           if (ir == n)  then
!              these are actually extrapolation points.
           ierr = ierr + next(2)
           go to 42
   41       continue
!           else
!              we should never have gotten here.
           go to 5005
!           end if
!        end if
   42    continue
!
     if (next(1) == 0)  go to 49
!        if (next(1) > 0)  then
!           in the current set of xe-points, there are next(1) to the
!           left of x(ir-1).
!
        if (ir > 2)  go to 43
!           if (ir == 2)  then
!              these are actually extrapolation points.
           ierr = ierr + next(1)
           go to 49
   43       continue
!           else
!              xe is not ordered relative to x, so must adjust
!              evaluation interval.
!
!              first, locate first point to left of x(ir-1).
!
           do i = jfirst, j-1
              if (xe(i) < x(ir-1))  go to 45
           end do
!
!              note-- cannot drop through here unless there is an error
!                     in chfdv.
           go to 5005
!
   45          continue
!              reset j.  (this will be the new jfirst.)
           j = i
!
!  now find out how far to back up in the x-array.
!
           do i = 1, ir-1
              if (xe(j) < x(i)) go to 47
           end do
!
!              nb-- can never drop through here, since xe(j)<x(ir-1).
!
   47          continue
!              at this point, either  xe(j) < x(1)
!                 or      x(i-1) <= xe(j) < x(i) .
!              reset ir, recognizing that it will be incremented before
!              cycling.
           ir = max(1, i-1)
!           end if
!        end if
   49    continue
!
     jfirst = j
!
!   end of ir-loop.
!
   50 continue
  ir = ir + 1
  if (ir <= n)  go to 10
!
!  normal return.
!
 5000 continue
  return
!
!  error returns.
!
 5001 continue
!     n<2 return.
  ierr = -1
  call xerror ('pchfd -- number of data points less than two', 44, ierr, 1)
  return

 5002 continue
!     incfd<1 return.
  ierr = -2
  call xerror ('pchfd -- increment less than one', 32, ierr, 1)
  return
!
 5003 continue
!     x-array not strictly increasing.
  ierr = -3
  call xerror ('pchfd -- x-array not strictly increasing', 40, ierr, 1)
  return
!
 5004 continue
!     ne<1 return.
  ierr = -4
  call xerror ('pchfd -- number of evaluation points less than one', 50, &
    ierr, 1)
  return
!
 5005 continue
!     error return from chfdv.
!  this case should never occur.
  ierr = -5
  call xerror ('pchfd -- error return from chfdv -- fatal', 41, ierr, 2)
  return
end
subroutine pchfe ( n, x, f, d, incfd, skip, ne, xe, fe, ierr )
!
!*******************************************************************************
!
!! PCHFE evaluates a piecewise cubic hermite function at an array of points.
!
!  may be used by itself for hermite interpolation,
!            or as an evaluator for pchim or pchic.
!
!  Description:
!
!          pchfe:  piecewise cubic hermite function evaluator
!
!     evaluates the cubic hermite function defined by  n, x, f, d  at
!     the points  xe(j), j=1(1)ne.
!
!     to provide compatibility with pchim and pchic, includes an
!     increment between successive values of the f- and d-arrays.
!
!  Parameters:
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of function values.  f(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     d -- (input) real array of derivative values.  d(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     incfd -- (input) increment between successive values in f and d.
!           (error return if  incfd<1 .)
!
!     skip -- (input/output) logical variable which should be set to
!           .true. if the user wishes to skip checks for validity of
!           preceding parameters, or to .false. otherwise.
!           this will save time in case these checks have already
!           been performed (say, in pchim or pchic).
!           skip will be set to .true. on normal return.
!
!     ne -- (input) number of evaluation points.  (error return if
!           ne<1 .)
!
!     xe -- (input) real array of points at which the function is to be
!           evaluated.
!
!          notes:
!           1. the evaluation will be most efficient if the elements
!              of xe are increasing relative to x;
!              that is,   xe(j) >= x(i)
!              implies    xe(k) >= x(i),  all k>=j .
!           2. if any of the xe are outside the interval [x(1),x(n)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     fe -- (output) real array of values of the cubic hermite function
!           defined by  n, x, f, d  at the points  xe.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning error:
!              ierr>0  means that extrapolation was performed at
!                 ierr points.
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!              ierr = -4  if ne<1 .
!             (the fe-array has not been changed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
!     2. most of the coding between the call to chfev and the end of
!        the ir-loop could be eliminated if it were permissible to
!        assume that xe is ordered relative to x.
!
!     3. chfev does not assume that x1 is less than x2.  thus, it would
!        be possible to write a version of pchfe that assumes a strict-
!        ly decreasing x-array by simply running the ir-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. the present code has a minor bug, which i have decided is not
!        worth the effort that would be required to fix it.
!        if xe contains points in [x(n-1),x(n)], followed by points <
!        x(n-1), followed by points >x(n), the extrapolation points
!        will be counted (at least) twice in the total returned in ierr.
!
  integer  n, incfd, ne, ierr
  real  x(n), f(incfd,n), d(incfd,n), xe(ne), fe(ne)
  logical  skip
  integer  i, ierc, ir, j, jfirst, next(2), nj
!
!  validity-check arguments.
!
  if (skip)  go to 5
!
  if ( n<2 )  go to 5001
  if ( incfd<1 )  go to 5002

  do i = 2, n
     if ( x(i)<=x(i-1) )  go to 5003
  end do
!
!function definition is ok, go on.
!
    5 continue
  if ( ne<1 )  go to 5004
  ierr = 0
  skip = .true.
!
!  loop over intervals.        (   interval index is  il = ir-1  . )
!                              ( interval is x(il)<=x<x(ir) . )
  jfirst = 1
  ir = 2
   10 continue
!
! skip out of loop if have processed all evaluation points.
!
     if (jfirst > ne)  go to 5000
!
!  locate all points in interval.
!
     do j = jfirst, ne
       if (xe(j) >= x(ir))  go to 30
     end do

     j = ne + 1
     go to 40
!
!  have located first point beyond interval.
!
   30    continue
     if (ir == n)  j = ne + 1

   40    continue
     nj = j - jfirst
!
!  skip evaluation if no points in interval.
!
     if (nj == 0)  go to 50
!
!  evaluate cubic at xe(i),  i = jfirst (1) j-1 .
!
    call chfev (x(ir-1),x(ir), f(1,ir-1),f(1,ir), d(1,ir-1),d(1,ir), &
      nj, xe(jfirst), fe(jfirst), next, ierc)

     if (ierc < 0)  go to 5005
!
     if (next(2) == 0)  go to 42
!        if (next(2) > 0)  then
!           in the current set of xe-points, there are next(2) to the
!           right of x(ir).
!
        if (ir < n)  go to 41
!           if (ir == n)  then
!              these are actually extrapolation points.
           ierr = ierr + next(2)
           go to 42
   41       continue
!           else
!              we should never have gotten here.
           go to 5005
!           end if
!        end if
   42    continue
!
     if (next(1) == 0)  go to 49
!        if (next(1) > 0)  then
!           in the current set of xe-points, there are next(1) to the
!           left of x(ir-1).
!
        if (ir > 2)  go to 43
!           if (ir == 2)  then
!              these are actually extrapolation points.
           ierr = ierr + next(1)
           go to 49
   43       continue
!           else
!              xe is not ordered relative to x, so must adjust
!              evaluation interval.
!
!              first, locate first point to left of x(ir-1).
           do i = jfirst, j-1
              if (xe(i) < x(ir-1))  go to 45
           end do
!              note-- cannot drop through here unless there is an error
!                     in chfev.
           go to 5005
!
   45          continue
!              reset j.  (this will be the new jfirst.)
           j = i
!
!              now find out how far to back up in the x-array.
           do i = 1, ir-1
              if (xe(j) < x(i)) go to 47
           end do
!              nb-- can never drop through here, since xe(j)<x(ir-1).
!
   47          continue
!              at this point, either  xe(j) < x(1)
!                 or      x(i-1) <= xe(j) < x(i) .
!              reset ir, recognizing that it will be incremented before
!              cycling.
           ir = max(1, i-1)
!           end if
!        end if
   49    continue
!
     jfirst = j
!
!   end of ir-loop.
!
   50 continue
  ir = ir + 1
  if (ir <= n)  go to 10
!
!  normal return.
!
 5000 continue
  return
!
!  error returns.
!
 5001 continue
!     n<2 return.
  ierr = -1
  call xerror ('pchfe -- number of data points less than two', 44, ierr, 1)
  return
!
 5002 continue
!     incfd<1 return.
  ierr = -2
  call xerror ('pchfe -- increment less than one', 32, ierr, 1)
  return
!
 5003 continue
!     x-array not strictly increasing.
  ierr = -3
  call xerror ('pchfe -- x-array not strictly increasing', 40, ierr, 1)
  return
!
 5004 continue
!     ne<1 return.
  ierr = -4
  call xerror ('pchfe -- number of evaluation points less than one', 50, ierr, 1)
  return
!
 5005 continue
!     error return from chfev.
!   *** this case should never occur.
  ierr = -5
  call xerror ('pchfe -- error return from chfev -- fatal', 41, ierr, 2)
  return
end
function pchia(n,x,f,d,incfd,skip,a,b,ierr)
!
!*******************************************************************************
!
!! PCHIA evaluates the integral of a piecewise cubic Hermite function.
!
!
!  Description:
!
!          pchia:  piecewise cubic hermite integrator, arbitrary limits
!
!     evaluates the definite integral of the cubic hermite function
!     defined by  n, x, f, d  over the interval [a, b].
!
!     to provide compatibility with pchim and pchic, includes an
!     increment between successive values of the f- and d-arrays.
!
!  Parameters:
!
!     value -- (output) value of the requested integral.
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of function values.  f(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     d -- (input) real array of derivative values.  d(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     incfd -- (input) increment between successive values in f and d.
!           (error return if  incfd<1 .)
!
!     skip -- (input/output) logical variable which should be set to
!           .true. if the user wishes to skip checks for validity of
!           preceding parameters, or to .false. otherwise.
!           this will save time in case these checks have already
!           been performed (say, in pchim or pchic).
!           skip will be set to .true. on return with ierr>=0 .
!
!     a,b -- (input) the limits of integration.
!           note:  there is no requirement that [a,b] be contained in
!                  [x(1),x(n)].  however, the resulting integral value
!                  will be highly suspect, if not.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning errors:
!              ierr = 1  if  a  is outside the interval [x(1),x(n)].
!              ierr = 2  if  b  is outside the interval [x(1),x(n)].
!              ierr = 3  if both of the above are true.  (note that this
!                        means that either [a,b] contains data interval
!                        or the intervals do not intersect at all.)
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!                (value has not been computed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
  integer incfd
  integer n
!
  real a
  real b
  real chfiv
  real d(incfd,n)
  real f(incfd,n)
  integer i
  integer ia
  integer ib
  integer ierd
  integer ierr
  integer ierv
  integer il
  integer ir
  real pchia
  real pchid
  logical skip
  real value
  real x(n)
  real xa
  real xb
!
!  validity-check arguments.
!
  if (skip)  go to 5
!
  if ( n<2 )  go to 5001
  if ( incfd<1 )  go to 5002

  do i = 2, n
    if ( x(i)<=x(i-1) )  go to 5003
  end do
!
!function definition is ok, go on.
!
    5 continue
  skip = .true.
  ierr = 0
  if ( (a<x(1)) .or. (a>x(n)) )  ierr = ierr + 1
  if ( (b<x(1)) .or. (b>x(n)) )  ierr = ierr + 2
!
!  compute integral value.
!
  if (a == b)  then
     value = 0.0E+00
  else
     xa = min (a, b)
     xb = max (a, b)
     if (xb <= x(2))  then
!           interval is to left of x(2), so use first cubic.
!
        value = chfiv ( x(1), x(2), f(1,1),f(1,2), &
          d(1,1),d(1,2), a, b, ierv)

        if (ierv < 0)  go to 5004
     else if (xa >= x(n-1))  then
!           interval is to right of x(n-1), so use last cubic.
!
        value = chfiv(x(n-1),x(n), f(1,n-1),f(1,n), &
          d(1,n-1),d(1,n), a, b, ierv)

        if (ierv < 0)  go to 5004
     else
!           'normal' case -- xa<xb, xa<x(n-1), xb>x(2).
!      ......locate ia and ib such that
!               x(ia-1)<xa<=x(ia)<=x(ib)<=xb<=x(ib+1)
        ia = 1
        do i = 1, n-1
           if (xa > x(i))  ia = i + 1
        end do
!             ia = 1 implies xa<x(1) .  otherwise,
!             ia is largest index such that x(ia-1)<xa,.
!
        ib = n
        do i = n, ia, -1
           if (xb < x(i))  ib = i - 1
        end do
!             ib = n implies xb>x(n) .  otherwise,
!             ib is smallest index such that xb<x(ib+1) .
!
!  Compute the integral.
!
        ierv = 0
        if (ib < ia)  then
!              this means ib = ia-1 and
!                 (a,b) is a subset of (x(ib),x(ia)).
!
           value = chfiv (x(ib),x(ia), f(1,ib),f(1,ia), &
             d(1,ib),d(1,ia), a, b, ierv)

           if (ierv < 0)  go to 5004
        else
!
!              first compute integral over (x(ia),x(ib)).
           if (ib == ia)  then
              value = 0.0E+00
           else

              value = pchid (n, x, f, d, incfd, skip, ia, ib, ierd)
!
              if (ierd < 0)  go to 5005
           end if
!
!              then add on integral over (xa,x(ia)).
           if (xa < x(ia))  then
              il = max (1, ia-1)
              ir = il + 1

              value = value + chfiv (x(il),x(ir), f(1,il),f(1,ir), &
                d(1,il),d(1,ir), xa, x(ia), ierv)

              if (ierv < 0)  go to 5004
           end if
!
!              then add on integral over (x(ib),xb).
           if (xb > x(ib))  then
              ir = min (ib+1, n)
              il = ir - 1

              value = value + chfiv (x(il),x(ir), f(1,il),f(1,ir), &
                d(1,il),d(1,ir), x(ib), xb, ierv)

              if (ierv < 0)  go to 5004
           end if
!
!              finally, adjust sign if necessary.
           if (a > b)  value = -value
        end if
     end if
  end if
!
!  normal return.
!
  pchia = value
  return
!
!  error returns.
!
 5001 continue
!     n<2 return.
  ierr = -1
  call xerror ('pchia -- number of data points less than two', 44, ierr, 1)
  return
!
 5002 continue
!     incfd<1 return.
  ierr = -2
  call xerror ('pchia -- increment less than one', 32, ierr, 1)
  return
!
 5003 continue
!     x-array not strictly increasing.
  ierr = -3
  call xerror ('pchia -- x-array not strictly increasing', 40, ierr, 1)
  return
!
 5004 continue
!     trouble in chfiv.  (should never occur.)
  ierr = -4
  call xerror ('pchia -- trouble in chfiv', 25, ierr, 1)
  return
!
 5005 continue
!     trouble in pchid.  (should never occur.)
  ierr = -5
  call xerror ('pchia -- trouble in pchid', 25, ierr, 1)
  return
end
subroutine pchic(ic,vc,switch,n,x,f,d,incfd,wk,nwk,ierr)
!
!*******************************************************************************
!
!! PCHIC sets derivatives needed to determine a piecewise monotone
!            piecewise cubic hermite interpolant to given data.
!
!            user control is available over boundary conditions and/or
!            treatment of points where monotonicity switches direction.
!
!  Description:
!
!         pchic:  piecewise cubic hermite interpolation coefficients.
!
!     sets derivatives needed to determine a piecewise monotone piece-
!     wise cubic interpolant to the data given in x and f satisfying the
!     boundary conditions specified by ic and vc.
!
!     the treatment of points where monotonicity switches direction is
!     controlled by argument switch.
!
!     to facilitate two-dimensional applications, includes an increment
!     between successive values of the f- and d-arrays.
!
!     the resulting piecewise cubic hermite function may be evaluated
!     by pchfe or pchfd.
!
!  Parameters:
!
!     ic -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           ic(1) = ibeg, desired condition at beginning of data.
!           ic(2) = iend, desired condition at end of data.
!
!           ibeg = 0  for the default boundary condition (the same as
!                     used by pchim).
!           if ibeg/=0, then its sign indicates whether the boundary
!                     derivative is to be adjusted, if necessary, to be
!                     compatible with monotonicity:
!              ibeg>0  if no adjustment is to be performed.
!              ibeg<0  if the derivative is to be adjusted for
!                     monotonicity.
!
!           allowable values for the magnitude of ibeg are:
!           ibeg = 1  if first derivative at x(1) is given in vc(1).
!           ibeg = 2  if second derivative at x(1) is given in vc(1).
!           ibeg = 3  to use the 3-point difference formula for d(1).
!                     (reverts to the default b.c. if n<3 .)
!           ibeg = 4  to use the 4-point difference formula for d(1).
!                     (reverts to the default b.c. if n<4 .)
!           ibeg = 5  to set d(1) so that the second derivative is con-
!              tinuous at x(2). (reverts to the default b.c. if n<4.)
!              this option is somewhat analogous to the "not a knot"
!              boundary condition provided by pchsp.
!
!          notes (ibeg):
!           1. an error return is taken if abs(ibeg)>5 .
!           2. only in case  ibeg<=0  is it guaranteed that the
!              interpolant will be monotonic in the first interval.
!              if the returned value of d(1) lies between zero and
!              3*slope(1), the interpolant will be monotonic.  this
!              is **not** checked if ibeg>0 .
!           3. if ibeg<0 and d(1) had to be changed to achieve mono-
!              tonicity, a warning error is returned.
!
!           iend may take on the same values as ibeg, but applied to
!           derivative at x(n).  in case iend = 1 or 2, the value is
!           given in vc(2).
!
!          notes (iend):
!           1. an error return is taken if abs(iend)>5 .
!           2. only in case  iend<=0  is it guaranteed that the
!              interpolant will be monotonic in the last interval.
!              if the returned value of d(1+(n-1)*incfd) lies between
!              zero and 3*slope(n-1), the interpolant will be monotonic.
!              this is **not** checked if iend>0 .
!           3. if iend<0 and d(1+(n-1)*incfd) had to be changed to
!              achieve monotonicity, a warning error is returned.
!
!     vc -- (input) real array of length 2 specifying desired boundary
!           values, as indicated above.
!           vc(1) need be set only if ic(1) = 1 or 2 .
!           vc(2) need be set only if ic(2) = 1 or 2 .
!
!     switch -- (input) indicates desired treatment of points where
!           direction of monotonicity switches:
!           set switch to zero if interpolant is required to be mono-
!           tonic in each interval, regardless of monotonicity of data.
!             notes:
!              1. this will cause d to be set to zero at all switch
!                 points, thus forcing extrema there.
!              2. the result of using this option with the default boun-
!                 dary conditions will be identical to using pchim, but
!                 will generally cost more compute time.
!                 this option is provided only to facilitate comparison
!                 of different switch and/or boundary conditions.
!           set switch nonzero to use a formula based on the 3-point
!              difference formula in the vicinity of switch points.
!           if switch is positive, the interpolant on each interval
!              containing an extremum is controlled to not deviate from
!              the data by more than switch*dfloc, where dfloc is the
!              maximum of the change of f on this interval and its two
!              immediate neighbors.
!           if switch is negative, no such control is to be imposed.
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of dependent variable values to be inter-
!           polated.  f(1+(i-1)*incfd) is value corresponding to x(i).
!
!     d -- (output) real array of derivative values at the data points.
!           these values will determine a monotone cubic hermite func-
!           tion on each subinterval on which the data are monotonic,
!           except possibly adjacent to switches in monotonicity.
!           the value corresponding to x(i) is stored in
!                d(1+(i-1)*incfd),  i=1(1)n.
!           no other entries in d are changed.
!
!     incfd -- (input) increment between successive values in f and d.
!           this argument is provided primarily for 2-d applications.
!           (error return if  incfd<1 .)
!
!     wk -- (scratch) real array of working storage.  the user may wish
!           to know that the returned values are:
!              wk(i)     = h(i)     = x(i+1) - x(i) ;
!              wk(n-1+i) = slope(i) = (f(1,i+1) - f(1,i)) / h(i)
!           for  i = 1(1)n-1.
!
!     nwk -- (input) length of work array.
!           (error return if  nwk<2*(n-1) .)
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning errors:
!              ierr = 1  if ibeg<0 and d(1) had to be adjusted for
!                        monotonicity.
!              ierr = 2  if iend<0 and d(1+(n-1)*incfd) had to be
!                        adjusted for monotonicity.
!              ierr = 3  if both of the above are true.
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!              ierr = -4  if abs(ibeg)>5 .
!              ierr = -5  if abs(iend)>5 .
!              ierr = -6  if both of the above are true.
!              ierr = -7  if nwk<2*(n-1) .
!             (the d-array has not been changed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
!  references  1. f.n.fritsch and r.e.carlson, 'monotone piecewise
!                 cubic interpolation,' siam j.numer.anal. 17, 2 (april
!                 1980), 238-246.
!               2. f.n.fritsch and j.butland, 'a method for constructing
!                 local monotone piecewise cubic interpolants,' llnl
!                 preprint ucrl-87559 (april 1982).
!               3. f.n.fritsch, 'piecewise cubic interpolation package,'
!                 llnl preprint ucrl-87285 (july 1982).
!
  integer incfd
  integer n
  integer nwk
!
  integer ic(2)
  integer ierr
  real switch
  real vc(2)
  real x(n), f(incfd,n), d(incfd,n), wk(nwk)
  integer  i, ibeg, iend, nless1
!
  if ( n<2 )  go to 5001
  if ( incfd<1 )  go to 5002

  do i = 2, n
    if ( x(i)<=x(i-1) )  go to 5003
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0
  if (abs(ibeg) > 5)  ierr = ierr - 1
  if (abs(iend) > 5)  ierr = ierr - 2
  if (ierr < 0)  go to 5004
!
!  function definition is ok -- go on.
!
  nless1 = n - 1
  if ( nwk < 2*nless1 )  go to 5007
!
!  set up h and slope arrays.
!
  do i = 1, nless1
    wk(i) = x(i+1) - x(i)
    wk(nless1+i) = (f(1,i+1) - f(1,i)) / wk(i)
  end do
!
!  special case n=2 -- use linear interpolation.
!
  if (nless1 > 1)  go to 1000
  d(1,1) = wk(2)
  d(1,n) = wk(2)
  go to 3000
!
!  normal case  (n >= 3) .
!
 1000 continue
!
!  set interior derivatives and default end conditions.
!
  call pchci (n, wk(1), wk(n), d, incfd)
!
!  set derivatives at points where monotonicity switches direction.
!
  if (switch == 0.0)  go to 3000

  call pchcs (switch, n, wk(1), wk(n), d, incfd, ierr)

  if (ierr /= 0)  go to 5008
!
!  set end conditions.
!
 3000 continue
  if ( (ibeg==0) .and. (iend==0) )  go to 5000

  call pchce (ic, vc, n, x, wk(1), wk(n), d, incfd, ierr)

  if (ierr < 0)  go to 5009
!
!  normal return.
!
 5000 continue
  return
!
!  error returns.
!
 5001 continue
!     n<2 return.
  ierr = -1
  call xerror ('pchic -- number of data points less than two', 44, ierr, 1)
  return

 5002 continue
!     incfd<1 return.
  ierr = -2
  call xerror ('pchic -- increment less than one', 32, ierr, 1)
  return

 5003 continue
!     x-array not strictly increasing.
  ierr = -3
  call xerror ('pchic -- x-array not strictly increasing', 40, ierr, 1)
  return
!
 5004 continue
!     ic out of range return.
  ierr = ierr - 3
  call xerror ('pchic -- ic out of range', 24, ierr, 1)
  return
!
 5007 continue
!     nwk < 2*(n-1)  return.
  ierr = -7
  call xerror ('pchic -- work array too small', 29, ierr, 1)
  return
!
 5008 continue
!     error return from pchcs.
  ierr = -8
  call xerror ('pchic -- error return from pchcs', 32, ierr, 1)
  return
!
 5009 continue
!     error return from pchce.
!   *** this case should never occur.
  ierr = -9
  call xerror ('pchic -- error return from pchce', 32, ierr, 1)
  return
end
function pchid(n,x,f,d,incfd,skip,ia,ib,ierr)
!
!*******************************************************************************
!
!! PCHID evaluates the definite integral of a piecewise cubic
!            hermite function over an interval whose endpoints are
!            data points.
!
!  Description:
!
!          pchid:  piecewise cubic hermite integrator, data limits
!
!     evaluates the definite integral of the cubic hermite function
!     defined by  n, x, f, d  over the interval [x(ia), x(ib)].
!
!     to provide compatibility with pchim and pchic, includes an
!     increment between successive values of the f- and d-arrays.
!
!  Parameters:
!
!     value -- (output) value of the requested integral.
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of function values.  f(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     d -- (input) real array of derivative values.  d(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     incfd -- (input) increment between successive values in f and d.
!           (error return if  incfd<1 .)
!
!     skip -- (input/output) logical variable which should be set to
!           .true. if the user wishes to skip checks for validity of
!           preceding parameters, or to .false. otherwise.
!           this will save time in case these checks have already
!           been performed (say, in pchim or pchic).
!           skip will be set to .true. on return with ierr = 0 or -4.
!
!     ia,ib -- (input) indices in x-array for the limits of integration.
!           both must be in the range [1,n].  (error return if not.)
!           no restrictions on their relative values.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!              ierr = -4  if ia or ib is out of range.
!                (value has not been computed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
  integer incfd
  integer  n
!
  real d(incfd,n)
  real f(incfd,n)
  real h
  integer i
  integer ia
  integer ib
  integer ierr
  integer iup
  integer low
  real pchid
  logical skip
  real sum
  real value
  real x(n)
!
!  validity-check arguments.
!
  if (skip)  go to 5
!
  if ( n<2 )  go to 5001
  if ( incfd<1 )  go to 5002

  do i = 2, n
    if ( x(i)<=x(i-1) )  go to 5003
  end do
!
!function definition is ok, go on.
!
    5 continue
  skip = .true.
  if ((ia<1) .or. (ia>n))  go to 5004
  if ((ib<1) .or. (ib>n))  go to 5004
  ierr = 0
!
!  compute integral value.
!
  if (ia == ib)  then
     value = 0.0E+00
  else
     low = min(ia, ib)
     iup = max(ia, ib) - 1
     sum = 0.0E+00
     do i = low, iup
        h = x(i+1) - x(i)
        sum = sum + h*( (f(1,i) + f(1,i+1)) + (d(1,i) - d(1,i+1))*(h/6.0) )
     end do
     value = 0.5 * sum
     if (ia > ib)  value = -value
  end if
!
!  normal return.
!
  pchid = value
  return
!
!  error returns.
!
 5001 continue
!     n<2 return.
  ierr = -1
  call xerror ('pchid -- number of data points less than two', 44, ierr, 1)
  return
!
 5002 continue
!     incfd<1 return.
  ierr = -2
  call xerror ('pchid -- increment less than one', 32, ierr, 1)
  return
!
 5003 continue
!     x-array not strictly increasing.
  ierr = -3
  call xerror ('pchid -- x-array not strictly increasing', 40, ierr, 1)
  return
!
 5004 continue
!     ia or ib out of range return.
  ierr = -4
  call xerror ('pchid -- ia or ib out of range', 30, ierr, 1)
  return
end
subroutine pchim ( n, x, f, d, incfd, ierr )
!
!*******************************************************************************
!
!! PCHIM sets derivatives for a piecewise cubic Hermite interpolant.
!
!
!  Discussion:
!
!    The routine set derivatives needed to determine a monotone piecewise
!    cubic hermite interpolant to given data.  the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.  (see pchic if user control is
!    desired over boundary or switch conditions.)
!
!    if the data are only piecewise monotonic, the interpolant will
!    have an extremum at each point where monotonicity switches direc-
!    tion.  (see pchic if user control is desired in such cases.)
!
!    to facilitate two-dimensional applications, includes an increment
!    between successive values of the f- and d-arrays.
!
!    the resulting piecewise cubic hermite function may be evaluated
!    by pchfe or pchfd.
!
!  References:
!
!    f.n.fritsch and r.e.carlson,
!    'monotone piecewise cubic interpolation,'
!    siam j.numer.anal.
!    17, 2 (april 1980), 238-246.
!
!    f.n.fritsch and j.butland,
!    'a method for constructing local monotone piecewise cubic interpolants,'
!    llnl preprint ucrl-87559 (april 1982).
!
!  Parameters:
!
!     n -- (input) number of data points.  (error return if n<2 .)
!           if n=2, simply does linear interpolation.
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of dependent variable values to be inter-
!           polated.  f(1+(i-1)*incfd) is value corresponding to x(i).
!           pchim is designed for monotonic data, but it will work for
!           any f-array.  it will force extrema at points where mono-
!           tonicity switches direction.  if some other treatment of
!           switch points is desired, pchic should be used instead.
!
!     d -- (output) real array of derivative values at the data points.
!           if the data are monotonic, these values will determine a
!           a monotone cubic hermite function.
!           the value corresponding to x(i) is stored in
!                d(1+(i-1)*incfd),  i=1(1)n.
!           no other entries in d are changed.
!
!     incfd -- (input) increment between successive values in f and d.
!           this argument is provided primarily for 2-d applications.
!           (error return if  incfd<1 .)
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning error:
!              ierr>0  means that ierr switches in the direction
!                 of monotonicity were detected.
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!             (the d-array has not been changed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
  integer incfd
  integer n
!
  real d(incfd,n)
  real del1
  real del2
  real dmax
  real dmin
  real drat1
  real drat2
  real dsave
  real f(incfd,n)
  real h1
  real h2
  real hsum
  real hsumt3
  integer i
  integer ierr
  integer nless1
  real pchst
  real w1
  real w2
  real x(n)
!
!  Check arguments.
!
  if ( n < 2 ) then
    ierr = -1
    call xerror ('pchim -- number of data points less than two', 44, ierr, 1)
    return
  end if

  if ( incfd < 1 ) then
    ierr = -2
    call xerror ('pchim -- increment less than one', 32, ierr, 1)
    return
  end if

  do i = 2, n
    if ( x(i)<=x(i-1) ) then
      ierr = -3
      call xerror ('pchim -- x-array not strictly increasing', 40, ierr, 1)
      return
    end if
  end do

  ierr = 0
  nless1 = n - 1
  h1 = x(2) - x(1)
  del1 = (f(1,2) - f(1,1))/h1
  dsave = del1
!
!  special case n=2 -- use linear interpolation.
!
  if ( n == 2 ) then
    d(1,1) = del1
    d(1,n) = del1
    return
  end if
!
!  normal case  (n >= 3).
!
  h2 = x(3) - x(2)
  del2 = (f(1,3) - f(1,2))/h2
!
!  set d(1) via non-centered three-point formula, adjusted to be
!  shape-preserving.
!
  hsum = h1 + h2
  w1 = (h1 + hsum)/hsum
  w2 = -h1/hsum
  d(1,1) = w1*del1 + w2*del2

  if ( pchst(d(1,1),del1) <= 0.0E+00 )  then
     d(1,1) = 0.0E+00
  else if ( pchst(del1,del2) < 0.0E+00 )  then
!        need do this check only if monotonicity switches.
     dmax = 3.0*del1
     if (abs(d(1,1)) > abs(dmax))  d(1,1) = dmax
  end if
!
!  loop through interior points.
!
  do i = 2, nless1

    if ( i > 2 ) then
      h1 = h2
      h2 = x(i+1) - x(i)
      hsum = h1 + h2
      del1 = del2
      del2 = (f(1,i+1) - f(1,i))/h2
    end if
!
!  set d(i)=0 unless data are strictly monotonic.
!
    d(1,i) = 0.0E+00
    if ( pchst(del1,del2) )  42, 41, 45
!
!  count number of changes in direction of monotonicity.
!
   41    continue
    if ( del2 /= 0.0E+00 ) then
      if ( pchst(dsave,del2) < 0.0E+00 )  ierr = ierr + 1
      dsave = del2
    end if
    go to 50

42  continue
    ierr = ierr + 1
    dsave = del2
    go to 50
!
!  use brodlie modification of butland formula.
!
45  continue

    hsumt3 = hsum+hsum+hsum
    w1 = (hsum + h1)/hsumt3
    w2 = (hsum + h2)/hsumt3
    dmax = max ( abs(del1), abs(del2) )
    dmin = min ( abs(del1), abs(del2) )
    drat1 = del1/dmax
    drat2 = del2/dmax
    d(1,i) = dmin/(w1*drat1 + w2*drat2)

   50 continue

  end do
!
!  set d(n) via non-centered three-point formula, adjusted to be
!  shape-preserving.
!
  w1 = -h2/hsum
  w2 = (h2 + hsum)/hsum
  d(1,n) = w1*del1 + w2*del2
  if ( pchst(d(1,n),del2) <= 0.0E+00 )  then
     d(1,n) = 0.0E+00
  else if ( pchst(del1,del2) < 0.0E+00 )  then
!        need do this check only if monotonicity switches.
     dmax = 3.0*del2
     if (abs(d(1,n)) > abs(dmax))  d(1,n) = dmax
  end if

  return
end
subroutine pchmc(n,x,f,d,incfd,skip,ismon,ierr)
!
!*******************************************************************************
!
!! PCHMC:  piecewise cubic hermite monotonicity checker.
!
!     checks the cubic hermite function defined by  n, x, f, d  for
!     monotonicity.
!
!     to provide compatibility with pchim and pchic, includes an
!     increment between successive values of the f- and d-arrays.
!
!  Parameters:
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of function values.  f(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     d -- (input) real array of derivative values.  d(1+(i-1)*incfd) is
!           the value corresponding to x(i).
!
!     incfd -- (input) increment between successive values in f and d.
!           (error return if  incfd<1 .)
!
!     skip -- (input/output) logical variable which should be set to
!           .true. if the user wishes to skip checks for validity of
!           preceding parameters, or to .false. otherwise.
!           this will save time in case these checks have already
!           been performed.
!           skip will be set to .true. on normal return.
!
!     ismon -- (output) integer array indicating on which intervals the
!           pch function defined by  n, x, f, d  is monotonic.
!           for data interval [x(i),x(i+1)],
!             ismon(i) = -1  if function is strictly decreasing;
!             ismon(i) =  0  if function is constant;
!             ismon(i) =  1  if function is strictly increasing;
!             ismon(i) =  2  if function is non-monotonic;
!             ismon(i) =  3  if unable to determine.  (this means that
!                            the d-values are near the boundary of the
!                            monotonicity region.  a small increase pro-
!                            duces non-monotonicity; decrease, strict
!                            monotonicity.)
!           the above applies to i=1(1)n-1.  ismon(n) indicates whether
!              the entire function is monotonic on [x(1),x(n)].
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!          (the ismon-array has not been changed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
!  references  f.n.fritsch and r.e.carlson, 'monotone piecewise cubic
!                 interpolation,' siam j.numer.anal. 17, 2 (april 1980),
!                 238-246.
!
  integer incfd
  integer  n
!
  real chfmc
  real d(incfd,n)
  real delta
  real f(incfd,n)
  integer i
  integer ierr
  integer ismon(n)
  integer nseg
  logical skip
  real x(n)
!
  if ( .not. skip ) then

    if ( n < 2 )  go to 5001
    if ( incfd<1 )  go to 5002

    do i = 2, n
      if ( x(i)<=x(i-1) )  go to 5003
    end do

    skip = .true.

  end if

  nseg = n - 1

  do i = 1, nseg

     delta = (f(1,i+1)-f(1,i))/(x(i+1)-x(i))

     ismon(i) = chfmc (d(1,i), d(1,i+1), delta)

     if (i == 1)  then
        ismon(n) = ismon(1)
     else
!           need to figure out cumulative monotonicity from following
!           'multiplication table'--
!
!                    *      i s m o n (i)
!                     *  -1   0   1   2   3
!               i      *--------------------*
!               s   -1 i -1  -1   2   2   3 i
!               m    0 i -1   0   1   2   3 i
!               o    1 i  2   1   1   2   3 i
!               n    2 i  2   2   2   2   2 i
!              (n)   3 i  3   3   3   2   3 i
!                      *--------------------*
!
!           if equal or already declared nonmonotonic, no change needed.
        if ((ismon(i)/=ismon(n)) .and. (ismon(n)/=2))  then
           if ( max(ismon(i), ismon(n)) > 1)  then
!                 at least one is either 'no' or 'maybe'.
              if (ismon(i) == 2)  then
                 ismon(n) = 2
              else
                 ismon(n) = 3
              end if
           else if (ismon(i)*ismon(n) < 0)  then
!                 both monotonic, but in opposite senses.
              ismon(n) = 2
           else
!                 at this point, one is zero, the other is +-1.
              ismon(n) = ismon(n) + ismon(i)
           end if
        end if
     end if

  end do

  ierr = 0
  return
!
!  error returns.
!
 5001 continue
!     n<2 return.
  ierr = -1
  call xerror ('pchmc -- number of data points less than two', 44, ierr, 1)
  return
!
 5002 continue
!     incfd<1 return.
  ierr = -2
  call xerror ('pchmc -- increment less than one', 32, ierr, 1)
  return
!
 5003 continue
!     x-array not strictly increasing.
  ierr = -3
  call xerror ('pchmc -- x-array not strictly increasing', 40, ierr, 1)
  return
end
function pchqa(n,x,f,d,a,b,ierr)
!
!*******************************************************************************
!
!! PCHQA: easy to use cubic hermite or spline integration
!             numerical integration, quadrature
!
!  purpose
!
!    evaluates the definite integral of a piecewise cubic hermite
!            or spline function over an arbitrary interval, easy to use.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!     evaluates the definite integral of the cubic hermite or spline
!     function defined by  n, x, f, d  over the interval [a, b].  this
!     is an easy to use driver for the routine pchia by f.n. fritsch
!     described in reference (2) below. that routine also has other
!     capabilities.
!
!  Parameters:
!
!     value -- (output) value of the requested integral.
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of function values.  f(i) is
!           the value corresponding to x(i).
!
!     d -- (input) real array of derivative values.  d(i) is
!           the value corresponding to x(i).
!
!     a,b -- (input) the limits of integration.
!           note:  there is no requirement that [a,b] be contained in
!                  [x(1),x(n)].  however, the resulting integral value
!                  will be highly suspect, if not.
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           warning errors:
!              ierr = 1  if  a  is outside the interval [x(1),x(n)].
!              ierr = 2  if  b  is outside the interval [x(1),x(n)].
!              ierr = 3  if both of the above are true.  (note that this
!                        means that either [a,b] contains data interval
!                        or the intervals do not intersect at all.)
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -3  if the x-array is not strictly increasing.
!                (value has not been computed in any of these cases.)
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!
!  references  1. f.n.fritsch and r.e.carlson, 'monotone piecewise
!                 cubic interpolation,' siam j.numer.anal. 17, 2 (april
!                 1980), 238-246.
!               2. f.n.fritsch, 'piecewise cubic hermite interpolation
!                 package, final specifications', lawrence livermore
!                 national laboratory, computer documentation ucid-30194,
!                 august 1982.
!
  integer n
!
  real a
  real b
  real d(n)
  real f(n)
  integer ierr
  integer, save :: incfd = 1
  real pchia
  real pchqa
  logical, save :: skip = .true.
  real x(n)
!
  pchqa  =  pchia( n, x, f, d, incfd, skip, a, b, ierr )

  return
end
subroutine pchsp(ic,vc,n,x,f,d,incfd,wk,nwk,ierr)
!
!*******************************************************************************
!
!! PCHSP: set derivatives needed to determine the hermite represen-
!            tation of the cubic spline interpolant to given data, with
!            specified boundary conditions.
!
!  Description:
!
!          pchsp:   piecewise cubic hermite spline
!
!     computes the hermite representation of the cubic spline inter-
!     polant to the data given in x and f satisfying the boundary
!     conditions specified by ic and vc.
!
!     to facilitate two-dimensional applications, includes an increment
!     between successive values of the f- and d-arrays.
!
!     the resulting piecewise cubic hermite function may be evaluated
!     by pchfe or pchfd.
!
!     note:  this is a modified version of c. de boor's cubic spline
!            routine cubspl.
!
!  Parameters:
!
!     ic -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           ic(1) = ibeg, desired condition at beginning of data.
!           ic(2) = iend, desired condition at end of data.
!
!           ibeg = 0  to set d(1) so that the third derivative is con-
!              tinuous at x(2).  this is the "not a knot" condition
!              provided by de boor's cubic spline routine cubspl.
!              < this is the default boundary condition. >
!           ibeg = 1  if first derivative at x(1) is given in vc(1).
!           ibeg = 2  if second derivative at x(1) is given in vc(1).
!           ibeg = 3  to use the 3-point difference formula for d(1).
!                     (reverts to the default b.c. if n<3 .)
!           ibeg = 4  to use the 4-point difference formula for d(1).
!                     (reverts to the default b.c. if n<4 .)
!          notes:
!           1. an error return is taken if ibeg is out of range.
!           2. for the "natural" boundary condition, use ibeg=2 and
!              vc(1)=0.
!
!           iend may take on the same values as ibeg, but applied to
!           derivative at x(n).  in case iend = 1 or 2, the value is
!           given in vc(2).
!
!          notes:
!           1. an error return is taken if iend is out of range.
!           2. for the "natural" boundary condition, use iend=2 and
!              vc(2)=0.
!
!     vc -- (input) real array of length 2 specifying desired boundary
!           values, as indicated above.
!           vc(1) need be set only if ic(1) = 1 or 2 .
!           vc(2) need be set only if ic(2) = 1 or 2 .
!
!     n -- (input) number of data points.  (error return if n<2 .)
!
!     x -- (input) real array of independent variable values.  the
!           elements of x must be strictly increasing:
!                x(i-1) < x(i),  i = 2(1)n.
!           (error return if not.)
!
!     f -- (input) real array of dependent variable values to be inter-
!           polated.  f(1+(i-1)*incfd) is value corresponding to x(i).
!
!     d -- (output) real array of derivative values at the data points.
!           these values will determine the cubic spline interpolant
!           with the requested boundary conditions.
!           the value corresponding to x(i) is stored in
!                d(1+(i-1)*incfd),  i=1(1)n.
!           no other entries in d are changed.
!
!     incfd -- (input) increment between successive values in f and d.
!           this argument is provided primarily for 2-d applications.
!           (error return if  incfd<1 .)
!
!     wk -- (scratch) real array of working storage.
!
!     nwk -- (input) length of work array.
!           (error return if nwk<2*n .)
!
!     ierr -- (output) error flag.
!           normal return:
!              ierr = 0  (no errors).
!           "recoverable" errors:
!              ierr = -1  if n<2 .
!              ierr = -2  if incfd<1 .
!              ierr = -3  if the x-array is not strictly increasing.
!              ierr = -4  if ibeg<0 or ibeg>4 .
!              ierr = -5  if iend<0 of iend>4 .
!              ierr = -6  if both of the above are true.
!              ierr = -7  if nwk is too small.
!               note:  the above errors are checked in the order listed,
!                   and following arguments have **not** been validated.
!             (the d-array has not been changed in any of these cases.)
!              ierr = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (the d-array may have been changed in this case.)
!             (             do **not** use it!                )
!
!  references  carl de boor, a practical guide to splines, springer-
!                 verlag (new york, 1978), pp. 53-59.
!
  integer incfd
  integer n
!
  real d(incfd,n)
  real f(incfd,n)
  real g
  integer ibeg
  integer ic(2)
  integer iend
  integer ierr
  integer index
  integer j
  integer nwk
  real pchdf
  real stemp(3)
  real vc(2)
  real wk(2,n)
  real x(n)
  real xtemp(4)
!
  if ( n<2 )  go to 5001
  if ( incfd<1 )  go to 5002

  do j = 2, n
     if ( x(j)<=x(j-1) )  go to 5003
  end do

  ibeg = ic(1)
  iend = ic(2)
  ierr = 0
  if ( (ibeg<0).or.(ibeg>4) )  ierr = ierr - 1
  if ( (iend<0).or.(iend>4) )  ierr = ierr - 2
  if ( ierr<0 )  go to 5004
!
!function definition is ok -- go on.
!
  if ( nwk < 2*n )  go to 5007
!
!  compute first differences of x sequence and store in wk(1,.). also,
!  compute first divided difference of data and store in wk(2,.).
!
  do j=2,n
     wk(1,j) = x(j) - x(j-1)
     wk(2,j) = (f(1,j) - f(1,j-1))/wk(1,j)
  end do
!
!  set to default boundary conditions if n is too small.
!
  if ( ibeg>n )  ibeg = 0
  if ( iend>n )  iend = 0
!
!  set up for boundary conditions.
!
  if ( (ibeg==1).or.(ibeg==2) )  then
     d(1,1) = vc(1)
  else if (ibeg > 2)  then
!
!        pick up first ibeg points, in reverse order.
!
     do j = 1, ibeg
        index = ibeg-j+1
        xtemp(j) = x(index)
        if (j < ibeg)  stemp(j) = wk(2,index)
     end do

     d(1,1) = pchdf (ibeg, xtemp, stemp, ierr)
     if (ierr /= 0)  go to 5009
     ibeg = 1
  end if

  if ( (iend==1).or.(iend==2) )  then
     d(1,n) = vc(2)
  else if (iend > 2)  then
!
!        pick up last iend points.
!
     do j = 1, iend
        index = n-iend+j
        xtemp(j) = x(index)
        if (j < iend)  stemp(j) = wk(2,index+1)
     end do

     d(1,n) = pchdf (iend, xtemp, stemp, ierr)

     if (ierr /= 0)  go to 5009
     iend = 1
  end if
!
!  begin coding from cubspl
!
!  a tridiagonal linear system for the unknown slopes s(j) of
!  f  at x(j), j=1,...,n, is generated and then solved by gauss elim-
!  ination, with s(j) ending up in d(1,j), all j.
!     wk(1,.) and wk(2,.) are used for temporary storage.
!
!  construct first equation from first boundary condition, of the form
!             wk(2,1)*s(1) + wk(1,1)*s(2) = d(1,1)
!
  if (ibeg == 0)  then
     if (n == 2)  then
!           no condition at left end and n = 2.
        wk(2,1) = 1.0E+00
        wk(1,1) = 1.0E+00
        d(1,1) = 2.0E+00 * wk(2,2)
     else
!           not-a-knot condition at left end and n > 2.
        wk(2,1) = wk(1,3)
        wk(1,1) = wk(1,2) + wk(1,3)
        d(1,1) =((wk(1,2) + 2.0E+00 * wk(1,1))*wk(2,2)*wk(1,3) &
                             + wk(1,2)**2*wk(2,3)) / wk(1,1)
     end if
  else if (ibeg == 1)  then
!        slope prescribed at left end.
     wk(2,1) = 1.0E+00
     wk(1,1) = 0.0E+00
  else
!        second derivative prescribed at left end.
     wk(2,1) = 2.0E+00
     wk(1,1) = 1.0E+00
     d(1,1) = 3.0*wk(2,2) - 0.5*wk(1,2)*d(1,1)
  end if
!
!  if there are interior knots, generate the corresponding equations and
!  carry out the forward pass of gauss elimination, after which the j-th
!  equation reads    wk(2,j)*s(j) + wk(1,j)*s(j+1) = d(1,j).
!
  if (n-1 > 1)  then
    do j=2,n-1
        if (wk(2,j-1) == 0.0E+00 )  go to 5008
        g = -wk(1,j+1)/wk(2,j-1)
        d(1,j) = g*d(1,j-1) + 3.0*(wk(1,j)*wk(2,j+1) + wk(1,j+1)*wk(2,j))
        wk(2,j) = g*wk(1,j-1) + 2.0E+00 *(wk(1,j) + wk(1,j+1))
    end do
  end if
!
!  construct last equation from second boundary condition, of the form
!           (-g*wk(2,n-1))*s(n-1) + wk(2,n)*s(n) = d(1,n)
!
!     if slope is prescribed at right end, one can go directly to back-
!     substitution, since arrays happen to be set up just right for it
!     at this point.
  if (iend == 1)  go to 30
!
  if (iend == 0)  then
     if (n==2 .and. ibeg==0)  then
!           not-a-knot at right endpoint and at left endpoint and n = 2.
        d(1,2) = wk(2,2)
        go to 30
     else if ((n==2) .or. (n==3 .and. ibeg==0))  then
!           either (n=3 and not-a-knot also at left) or (n=2 and *not*
!           not-a-knot at left end point).
        d(1,n) = 2.0E+00 *wk(2,n)
        wk(2,n) = 1.0E+00
        if (wk(2,n-1) == 0.0E+00 )  go to 5008
        g = -1.0/wk(2,n-1)
     else
!           not-a-knot and n >= 3, and either n>3 or  also not-a-
!           knot at left end point.
        g = wk(1,n-1) + wk(1,n)
!           do not need to check following denominators (x-differences).
        d(1,n) = ((wk(1,n)+2.0E+00 *g)*wk(2,n)*wk(1,n-1) &
          + wk(1,n)**2*(f(1,n-1)-f(1,n-2))/wk(1,n-1))/g
        if (wk(2,n-1) == 0.0E+00 )  go to 5008
        g = -g/wk(2,n-1)
        wk(2,n) = wk(1,n-1)
     end if
  else
!        second derivative prescribed at right endpoint.
     d(1,n) = 3.0*wk(2,n) + 0.5*wk(1,n)*d(1,n)
     wk(2,n) = 2.0E+00
     if (wk(2,n-1) == 0.0E+00 )  go to 5008
     g = -1.0/wk(2,n-1)
  end if
!
!  complete forward pass of gauss elimination.
!
  wk(2,n) = g*wk(1,n-1) + wk(2,n)
  if (wk(2,n) == 0.0E+00 )   go to 5008
  d(1,n) = (g*d(1,n-1) + d(1,n))/wk(2,n)
!
!  carry out back substitution
!
   30 continue

  do j=n-1,1,-1
     if (wk(2,j) == 0.0E+00 )  go to 5008
     d(1,j) = (d(1,j) - wk(1,j)*d(1,j+1))/wk(2,j)
  end do

  return
!
!  error returns.
!
 5001 continue
!     n<2 return.
  ierr = -1
  call xerror ('pchsp -- number of data points less than two', 44, ierr, 1)
  return

 5002 continue
!     incfd<1 return.
  ierr = -2
  call xerror ('pchsp -- increment less than one', 32, ierr, 1)
  return

 5003 continue
!     x-array not strictly increasing.
  ierr = -3
  call xerror ('pchsp -- x-array not strictly increasing', 40, ierr, 1)
  return
!
 5004 continue
!     ic out of range return.
  ierr = ierr - 3
  call xerror ('pchsp -- ic out of range', 24, ierr, 1)
  return

 5007 continue
!     nwk too small return.
  ierr = -7
  call xerror ('pchsp -- work array too small', 29, ierr, 1)
  return
!
 5008 continue
!     singular system.
!   *** theoretically, this can only occur if successive x-values
!   *** are equal, which should already have been caught (ierr=-3).
  ierr = -8
  call xerror ('pchsp -- singular linear system', 31, ierr, 1)
  return
!
 5009 continue
!     error return from pchdf.
!   *** this case should never occur.
  ierr = -9
  call xerror ('pchsp -- error return from pchdf', 32, ierr, 1)
  return
end
function pchst(arg1,arg2)
!
!*******************************************************************************
!
!! PCHST: pchip sign-testing routine.
!
!
!     returns:
!        -1. if arg1 and arg2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if arg1 and arg2 are of the same sign.
!
!     the object is to do this without multiplying arg1*arg2, to avoid
!     possible over/underflow problems.
!
!  programmed by:  fred n. fritsch,  fts 532-4275, (415) 422-4275,
!                  mathematics and statistics division,
!                  lawrence livermore national laboratory.
!
!
  real arg1
  real arg2
  real pchst
!
!  perform the test.
!
  pchst = sign(1.0,arg1) * sign(1.0,arg2)
  if ((arg1==0.0) .or. (arg2==0.0))  pchst = 0.0E+00
!
  return
end
subroutine pchsw(dfmax,iextrm,d1,d2,h,slope,ierr)
!
!*******************************************************************************
!
!! PCHSW: pchcs switch excursion limiter.
!
!
!     called by  pchcs  to adjust d1 and d2 if necessary to insure that
!     the extremum on this interval is not further than dfmax from the
!     extreme data value.
!
!  Parameters:
!
!     dfmax -- (input) maximum allowed difference between f(iextrm) and
!           the cubic determined by derivative values d1,d2.  (assumes
!           dfmax>0.)
!
!     iextrm -- (input) index of the extreme data value.  (assumes
!           iextrm = 1 or 2 .  any value /=1 is treated as 2.)
!
!     d1,d2 -- (input) derivative values at the ends of the interval.
!           (assumes d1*d2 <= 0.)
!          (output) may be modified if necessary to meet the restriction
!           imposed by dfmax.
!
!     h -- (input) interval length.  (assumes  h>0.)
!
!     slope -- (input) data slope on the interval.
!
!     ierr -- (output) error flag.  should be zero.
!           if ierr=-1, assumption on d1 and d2 is not satisfied.
!           if ierr=-2, quadratic equation locating extremum has
!                       negative descriminant (should never occur).
!
  real cp
  real d1
  real d2
  real dfmax
  real dmax
  real, parameter :: fact = 100.0E+00
  real h
  integer ierr
  integer iextrm
  real lambda
  real nu
  real phi
  real radcal
  real rho
  real r1mach
  real sigma
  real slope
  real small
  real that
  real, parameter :: third = 0.33333
!
!  notation and general remarks.
!
!     rho is the ratio of the data slope to the derivative being tested.
!     lambda is the ratio of d2 to d1.
!     that = t-hat(rho) is the normalized location of the extremum.
!     phi is the normalized value of p(x)-f1 at x = xhat = x-hat(rho),
!           where  that = (xhat - x1)/h .
!        that is, p(xhat)-f1 = d*h*phi,  where d=d1 or d2.
!     similarly,  p(xhat)-f2 = d*h*(phi-rho) .
!
!      small should be a few orders of magnitude greater than macheps.
!
  small = fact * epsilon ( 1.0E+00 )
!
!  do main calculation.
!
  if (d1 == 0.0)  then
!
!        special case -- d1==0.0E+00 .
!
!          if d2 is also zero, this routine should not have been called.
     if (d2 == 0.0)  go to 5001
!
     rho = slope/d2
!          extremum is outside interval when rho >= 1/3 .
     if (rho >= third)  go to 5000
     that = (2.0*(3.0*rho-1.0)) / (3.0*(2.0*rho-1.0))
     phi = that**2 * ((3.0*rho-1.0)/3.0)
!
!          convert to distance from f2 if iextrm/=1 .
     if (iextrm /= 1)  phi = phi - rho
!
!          test for exceeding limit, and adjust accordingly.
     dmax = dfmax / (h*abs(phi))
     if (abs(d2) > dmax)  d2 = sign (dmax, d2)
  else
!
     rho = slope/d1
     lambda = -d2/d1
     if (d2 == 0.0)  then
!
!           special case -- d2==0.0E+00 .
!
!             extremum is outside interval when rho >= 1/3 .
        if (rho >= third)  go to 5000
        cp = 2.0E+00 - 3.0*rho
        nu = 1.0E+00 - 2.0*rho
        that = 1.0E+00 / (3.0*nu)
     else
        if (lambda <= 0.0)  go to 5001
!
!           normal case -- d1 and d2 both nonzero, opposite signs.
!
        nu = 1.0E+00 - lambda - 2.0*rho
        sigma = 1.0E+00 - rho
        cp = nu + sigma
        if (abs(nu) > small)  then
           radcal = (nu - (2.0*rho+1.0))*nu + sigma**2
           if (radcal < 0.0)  go to 5002
           that = (cp - sqrt(radcal)) / (3.0*nu)
        else
           that = 1.0/(2.0*sigma)
        end if
     end if
     phi = that*((nu*that - cp)*that + 1.0)
!
!          convert to distance from f2 if iextrm/=1 .
     if (iextrm /= 1)  phi = phi - rho
!
!          test for exceeding limit, and adjust accordingly.
     dmax = dfmax / (h*abs(phi))
     if (abs(d1) > dmax)  then
        d1 = sign (dmax, d1)
        d2 = -lambda*d1
     end if
  end if
!
!  normal return.
!
 5000 continue
  ierr = 0
  return
!
!  error returns.
!
 5001 continue
!     d1 and d2 both zero, or both nonzero and same sign.
  ierr = -1
  call xerror ('pchsw -- d1 and/or d2 invalid', 29, ierr, 1)
  return
!
 5002 continue
!     negative value of radical (should never occur).
  ierr = -2
  call xerror ('pchsw -- negative radical', 25, ierr, 1)

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
function pimach ()
!
!*******************************************************************************
!
!! PIMACH returns the value of pi.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Output, real PIMACH, the value of PI.
!
  real pimach
!
  pimach = 4.0E+00 * atan ( 1.0E+00 )

  return
end
subroutine q1da(f,a,b,eps,r,e,kf,iflag)
!
!*******************************************************************************
!
!! Q1DA approximates the definite integral of a user defined function of one variable.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!       f     (input) the name of the external function which
!             evaluates the integrand.
!       a
!       b     (input) the endpoints of the integration interval
!       eps   (input) the accuracy to which you want the integral
!                computed.  if you want 2 digits of accuracy set
!                eps=.01, for 3 digits set eps=.001, etc.
!                eps must be positive.
!       r     (output) q1da's best estimate of your integral
!       e     (output) an estimate of abs(integral-r)
!       kf    (output) the cost of the integration, measured in
!                   number of evaluations of your integrand.
!                   kf will always be at least 30.
!       iflag (output) termination flag...possible values are
!               0   normal completion, e satisfies
!                        e<eps  and  e<eps*abs(r)
!               1   normal completion, e satisfies
!                        e<eps, but e>eps*abs(r)
!               2   normal completion, e satisfies
!                        e<eps*abs(r), but e>eps
!               3   normal completion but eps was too small to
!                     satisfy absolute or relative error request.
!
!               4   aborted calculation because of serious rounding
!                     error.  probably e and r are consistent.
!               5   aborted calculation because of insufficient storage.
!                     r and e are consistent.
!               6   aborted calculation because of serious difficulties
!                     meeting your error request.
!               7   aborted calculation because eps was set <= 0.0E+00
!
!            note...if iflag=3, 4, 5 or 6 consider using q1dax instead.
!
!    w h e r e   i s   y o u r   i n t e g r a n d ?
!
!        you must write a fortran function, called f, to evaluate
!        the integrand.  usually this looks like...
!               function f(x)
!                    f=(evaluate the integrand at the point x)
!                    return
!               end
!
!
!    t y p i c a l   p r o b l e m   s e t u p
!
!          a=0.0E+00
!          b=1.0E+00          (set interval endpoints to [0,1])
!          eps=0.001       (set accuracy request for 3 digits)
!          call q1da(f,a,b,eps,r,e,kf,iflag)
!        end
!        function f(x)
!            f=sin(2.*x)-sqrt(x)     (for example)
!            return
!        end
!      for this sample problem, the output is
!  0.0E+00    1.0E+00     .001    .041406750    .69077e-07    30    0
!
!    r e m a r k   i.
!
!           a small amout of randomization is built into this program.
!           calling q1da a few times in succession will give different
!           but hopefully consistent results.
!
!   r e m a r k   ii.
!
!           this routine is designed for integration over a finite
!           interval.  thus the input arguments a and b must be
!           valid real numbers on your computer.  if you want to do
!           an integral over an infinite interval set a or b or both
!           large enough so that the interval [a,b] contains most of
!           the integrand.  care is necessary, however.  for example,
!           to integrate exp(-x*x) on the entire real line one could
!           take a=-20., b=20. or similar values to get good results.
!           if you took a=-1.e10 and b=+1.e10 two bad things would
!           occur. first, you will certainly get an error message from
!           the exp routine, as its argument is too small.  other
!           things could happen too, for example an underflow.
!           second, even if the arithmetic worked properly q1da will
!           surely give an incorrect answer, because its first try
!           at sampling the integrand is based on your scaling and
!           it is very unlikely to select evaluation points in the
!           infinitesmally small interval [-20,20] where all the
!           integrand is concentrated, when a, b are so large.
!
!    m o r e   f l e x i b i l i t y
!
!           q1da is an easy to use driver for another program, q1dax.
!           q1dax provides several options which are not available
!                with q1da.
!
  integer, parameter :: nmax = 50
!
  real a
  real b
  real e
  real eps
  real f
  real fmax
  real fmin
  integer iflag
  integer kf
  integer nint
  real r
  logical rst
  real w(nmax,6)
!
  external f
!
  nint=1
  rst = .false.
  call q1dax(f,a,b,eps,r,e,nint,rst,w,nmax,fmin,fmax,kf,iflag)

  return
end
subroutine q1dax(f,a,b,eps,r,e,nint,rst,w,nmax,fmin,fmax,kf,iflag)
!
!*******************************************************************************
!
!! Q1DAX approximates the integral of a function of one variable.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!             for an easier to use routine see q1da.
!
!             capabilities of q1dax (in addition to those of q1da)
!             include:
!                ability to restart a calculation to greater
!                  accuracy without penalty...
!                ability to specify an initial partition of
!                  the integration interval...
!                ability to increase the work space to handle
!                  more difficult problems...
!                output of largest/smallest integrand value for
!                  applications such as scaling graphs...
!
!       a r g u m e n t s   i n   t h e   c a l l   s e q u e n c e
!
!       f      (input) the name of your integrand function.
!                      this name must appear in an external statement
!                      in any program which calls q1dax.
!                      you must write f in the form
!                         function f(x)
!                              f=(evaluate integrand at the point x)
!                              return
!                         end
!       a
!       b       (input) endpoints of integration interval
!       eps     (input)  accuracy to which the integral is to be calculated.
!                          q1dax will try to achieve relative accuracy,
!                          e.g. set eps=.01 for 2 digits, .001 for 3, etc.
!       r       (output) the estimate of the integral
!       e       (output) the estimate of the absolute error in r.
!       nint    (input
!                output)
!                        as an input quantity, nint must be set to
!                          the number of subintervals in the initial
!                          partition of [a,b].  for most problems
!                          this is just 1, the interval [a,b] itself.
!                          nint must be less than nmax, see below.
!                          nint is useful if you would like to help
!                          q1dax locate a difficult spot on [a,b].
!                          in this regard nint is used along
!                          with the array w (see below).  if you set
!                          nint=1 it is not necessary to be concerned
!                          with w, except that it must be dimensioned...
!                          as an example of more general applications,
!                          if [a,b]=[0,1] but the integrand jumps at 0.3,
!                          it would be wise to set nint=2 and then set
!                                  w(1,1)=0.0E+00  (left endpoint)
!                                  w(2,1)=0.3  (singular point)
!                                  w(3,1)=1.0E+00  (right endpoint)
!                          if you set nint greater than 1, be sure to
!                          check that you have also set
!                                w(1,1)=a  and  w(nint+1,1)=b
!                        as an output quantity, nint gives the
!                          number of subintervals in the final
!                          partition of [a,b].
!       rst     (input) a logical variable (e.g. true or false)
!                       set rst=.false. for initial call to q1dax
!                       set rst=.true. for a subsequent call,
!                           e.g. one for which more accuracy is
!                           desired (smaller eps).  a restart only
!                           makes sense if the preceding call returned
!                           with a value of iflag (see below) less than 3.
!                           on a restart you may not change the values of any
!                           other arguments in the call sequence, except eps.
!       w(nmax,6) (input and
!                  scratch)
!                       w is an array used by q1dax.
!                        you  m u s t  include a dimension statement in
!                        your calling program to allocate this storage.
!                        this should be of the form
!                                   dimension w(nmax,6)
!                        where nmax is an integer. an adequate value of
!                        nmax is 50.  if you set nint>1 you must also
!                        initialize w, see nint above.
!       nmax    (input) an integer equal to the first subscript in the
!                        dimension statement for the array w.  this is
!                        also equal to the maximum number of subintervals
!                        permitted in the internal partition of [a,b].
!                        a value of 50 is ample for most problems.
!       fmin
!       fmax    (output) the smallest and largest values of the integrand
!                          which occurred during the calculation.  the
!                          actual integrand range on [a,b] may, of course,
!                          be greater but probably not by more than 10%.
!       kf      (output) the actual number of integrand evaluations used
!                          by q1dax to approximate this integral.  kf
!                          will always be at least 30.
!       iflag   (output) termination flag...possible values are
!                  0   normal completion, e satisfies
!                           e<eps  and  e<eps*abs(r)
!                  1   normal completion, e satisfies
!                           e<eps, but e>eps*abs(r)
!                  2   normal completion, e satisfies
!                           e<eps*abs(r), but e>eps
!                  3   normal completion but eps was too small to
!                        satisfy absolute or relative error request.
!                  4   aborted calculation because of serious rounding
!                        error.  probably e and r are consistent.
!                  5   aborted calculation because of insufficient storage.
!                        r and e are consistent.  perhaps increasing nmax
!                        will produce better results.
!                  6   aborted calculation because of serious difficulties
!                        meeting your error request.
!                  7   aborted calculation because either eps, nint or nmax
!                        has been set to an illegal value.
!                  8   aborted calculation because you set nint>1 but forgot
!                        to set w(1,1)=a  and  w(nint+1,1)=b
!
!     t y p i c a l   p r o b l e m   s e t   u p
!
!      dimension w(50,6)
!      logical rst
!      external f
!      a=0.0E+00
!      b=1.0E+00
!      w(1,1)=a
!      w(2,1)=.3      [set internal partition point at .3]
!      w(3,1)=b
!      nint=2         [initial partition has 2 intervals]
!      rst=.false.
!      eps=.001
!      nmax=50
!
!    1 call q1dax(f,a,b,eps,r,e,nint,rst,w,nmax,fmin,fmax,kf,iflag)
!
!      if ( eps== .0001 .or. iflag>=3)stop
!      rst=.true
!      eps=.0001      [ask for another digit]
!      go to 1
!    end
!    function f(x)
!      if ( x< .3)
!     1  then
!          f=x**(0.2)*log(x)
!        else
!          f=sin(x)
!      end if
!      return
!    end
!
!
!            r e m a r k
!
!               when you use q1adx with nint=1, we have built a small
!            amount of randomization into it.  repeated calls during
!            the same run will produce different, but hopefully
!            consistent, results.
!
  integer nmax
!
  real a
  real b
  integer c
  real e
  real eb
  real epmach
  real eps
  real f
  real fmax
  real fmaxl
  real fmaxr
  real fmin
  real fminl
  real fminr
  real fmn
  real fmx
  integer i
  integer iflag
  integer iroff
  integer isamax
  integer kf
  integer loc
  integer mxtry
  integer nint
  real r
  real r1mach
  real rab
  real rabs
  real rav
  logical rst
  real t
  real te
  real te1
  real te2
  real tr
  real tr1
  real tr2
  real uflow
  real uni
  real w(nmax,6)
  real xm
!
  external f
!
  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  mxtry=nmax/2
!
!  in case there is no more room, we can toss out easy intervals,
!  at most mxtry times.
!
  if ( a==b) then
    r=0.0E+00
    e=0.0E+00
    nint=0
    iflag=0
    kf=1
    fmin=f(a)
    fmax=fmin
    go to 20
  end if

  if ( rst) then
     if ( iflag<3) then
       eb=max(100.*uflow,max(eps,50.*epmach)*abs(r))
       do i=1,nint
           if ( abs(w(i,3))>(eb*(w(i,2)-w(i,1))/(b-a)) ) then
                             w(i,3)=abs(w(i,3))
           else
                             w(i,3)=-abs(w(i,3))
           end if
       end do
       go to 15
     else
       go to 20
     end if
  end if

  kf=0

  if ( eps <= 0. .or. nint <= 0 .or. nint >= nmax) then
      iflag=7
      go to 20
  end if

  if ( nint==1) then
      w(1,1)=a
      w(2,2)=b
      w(1,5)=a
      w(1,6)=b
      w(2,5)=a
      w(2,6)=b
!
!  Select the first subdivision randomly.
!
      w(1,2) = a + ( b - a ) / 2.0 * ( 2.0 * uni() + 7.0 ) / 8.0
      w(2,1) = w(1,2)
      nint = 2
    else
      if ( w(1,1)/=a .or. w(nint+1,1)/=b) then
           iflag=8
           go to 20
      end if
      w(1,5)=a
      do i=1,nint
         w(i,2)=w(i+1,1)
         w(i,5)=w(i,1)
         w(i,6)=w(i,2)
      end do
  end if

  iflag = 0
  iroff=0
  rabs=0.0E+00

  do i=1,nint

      call gl15t(f,w(i,1),w(i,2),dble(w(i,5)),dble(w(i,6)), &
               w(i,4),w(i,3),rab,rav,fmn,fmx)
      kf=kf+15

      if ( i==1) then
        r=w(i,4)
        e=w(i,3)
        rabs=rabs+rab
        fmin=fmn
        fmax=fmx
      else
        r=r+w(i,4)
        e=e+w(i,3)
        rabs=rabs+rab
        fmax=max(fmax,fmx)
        fmin=min(fmin,fmn)
      end if

  end do

  w(nint+1:nmax,3) = 0.0E+00

   15 continue
!
!   main subprogram loop
!
  if ( 100.*epmach*rabs>=abs(r) .and. e<eps)go to 20
  eb=max(100.*uflow,max(eps,50.*epmach)*abs(r))
  if ( e<=eb) go to 20

  if (nint<nmax) then
    nint = nint+1
    c = nint
  else
    c=0
16  continue
    if ( c==nmax .or. mxtry<=0) then
        iflag=5
        go to 20
    end if

    c=c+1
    if ( w(c,3)>0.0) go to 16
!
!  found an interval to throw out
!
    mxtry=mxtry-1
  end if

  loc=isamax(nint,w(1,3),1)
  xm = w(loc,1)+(w(loc,2)-w(loc,1))/2.
  if ((max(abs(w(loc,1)),abs(w(loc,2))))> &
    ((1.+100.*epmach)*(abs(xm)+0.1e+04*uflow))) then

        call gl15t(f,w(loc,1),xm,dble(w(loc,5)),dble(w(loc,6)), &
          tr1,te1,rab,rav,fminl,fmaxl)

        kf=kf+15
        if (te1<(eb*(xm-w(loc,1))/(b-a))) te1=-te1
        call gl15t(f,xm,w(loc,2),dble(w(loc,5)),dble(w(loc,6)), &
          tr2,te2,rab,rav,fminr,fmaxr)

        kf=kf+15
        fmin=min(fmin,fminl,fminr)
        fmax=max(fmax,fmaxl,fmaxr)
        if (te2<(eb*(w(loc,2)-xm)/(b-a))) te2=-te2
        te = abs(w(loc,3))
        tr = w(loc,4)
        w(c,3) = te2
        w(c,4) = tr2
        w(c,1) = xm
        w(c,2) = w(loc,2)
        w(c,5) = w(loc,5)
        w(c,6) = w(loc,6)
        w(loc,3) = te1
        w(loc,4) = tr1
        w(loc,2) = xm
        e = e-te+(abs(te1)+abs(te2))
        r = r-tr+(tr1+tr2)
        if ( abs(abs(te1)+abs(te2)-te)<0.001*te) then
            iroff=iroff+1
            if ( iroff>=10) then
                 iflag=4
                 go to 20
            end if
        end if
      else
        if (eb>w(loc,3))then
               w(loc,3) = 0.
           else
               iflag=6
               go to 20
        end if
  end if
  go to 15
!
!  all exits from here
!
   20 continue
  if ( iflag>=4)return
  iflag=3
  t=eps*abs(r)
  if ( e>eps .and. e>t)return
  iflag=2
  if ( e>eps .and. e<t)return
  iflag=1
  if ( e<eps .and. e>t)return
  iflag=0
  return
end
subroutine qagi(f,bound,inf,epsabs,epsrel,result,abserr,neval, &
  ier,limit,lenw,last,iwork,work)
!
!*******************************************************************************
!
!! QAGI calculates an approximation result to a given
!            integral   i = integral of f over (bound,+infinity)
!                    or i = integral of f over (-infinity,bound)
!                    or i = integral of f over (-infinity,+infinity)
!            hopefully satisfying following claim for accuracy
!            abs(i-result)<=max(epsabs,epsrel*abs(i)).
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!         on entry
!            f      - real
!                   function subprogram defining the integrand
!                   function f(x). the actual name for f needs to be
!                     declared e x t e r n a l in the driver program.
!
!            bound  - real
!                     finite bound of integration range
!                     (has no meaning if interval is doubly-infinite)
!
!            inf    - integer
!                     indicating the kind of integration range involved
!                     inf = 1 corresponds to  (bound,+infinity),
!                     inf = -1            to  (-infinity,bound),
!                     inf = 2             to (-infinity,+infinity).
!
!            epsabs - real
!                     absolute accuracy requested
!            epsrel - real
!                     relative accuracy requested
!                     if  epsabs<=0
!                     and epsrel<max(50*rel.mach.acc.,0.5d-28),
!                     the routine will end with ier = 6.
!
!
!         on return
!            result - real
!                     approximation to the integral
!
!            abserr - real
!                     estimate of the modulus of the absolute error,
!                     which should equal or exceed abs(i-result)
!
!            neval  - integer
!                     number of integrand evaluations
!
!            ier    - integer
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                   - ier>0 abnormal termination of the routine. the
!                             estimates for result and error are less
!                             reliable. it is assumed that the requested
!                             accuracy has not been achieved.
!            error messages
!                     ier = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more
!                             subdivisions by increasing the value of
!                             limit (and taking the according dimension
!                             adjustments into account). however, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. if
!                             the position of a local difficulty can be
!                             determined (e.g. singularity,
!                             discontinuity within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. if possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 the occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 the algorithm does not converge.
!                             roundoff error is detected in the
!                             extrapolation table.
!                             it is assumed that the requested tolerance
!                             cannot be achieved, and that the returned
!                             result is the best which can be obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             (epsabs<=0 and
!                              epsrel<max(50*rel.mach.acc.,0.5d-28))
!                              or limit<1 or leniw<limit*4.
!                             result, abserr, neval, last are set to
!                             zero. exept when limit or leniw is
!                             invalid, iwork(1), work(limit*2+1) and
!                             work(limit*3+1) are set to zero, work(1)
!                             is set to a and work(limit+1) to b.
!
!         dimensioning parameters
!            limit - integer
!                    dimensioning parameter for iwork
!                    limit determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (a,b), limit>=1. in many cases limit = 100 is ok.
!                    if limit<1, the routine will end with ier = 6.
!
!            lenw  - integer
!                    dimensioning parameter for work
!                    lenw must be at least limit*4.
!                    if lenw<limit*4, the routine will end
!                    with ier = 6.
!
!            last  - integer
!                    on return, last equals the number of subintervals
!                    produced in the subdivision process, which
!                    determines the number of significant elements
!                    actually in the work arrays.
!
!         work arrays
!            iwork - integer
!                    vector of dimension at least limit, the first
!                    k elements of which contain pointers
!                    to the error estimates over the subintervals,
!                    such that work(limit*3+iwork(1)),... ,
!                    work(limit*3+iwork(k)) form a decreasing
!                    sequence, with k = last if last<=(limit/2+2), and
!                    k = limit+1-last otherwise
!
!            work  - real
!                    vector of dimension at least lenw
!                    on return
!                    work(1), ..., work(last) contain the left
!                   end points of the subintervals in the
!                     partition of (a,b),
!                    work(limit+1), ..., work(limit+last) contain
!                     the right end points,
!                    work(limit*2+1), ...,work(limit*2+last) contain the
!                     integral approximations over the subintervals,
!                    work(limit*3+1), ..., work(limit*3)
!                     contain the error estimates.
!
  integer lenw
  integer limit
!
  real abserr
  real bound
  real epsabs
  real epsrel
  real f
  integer ier
  integer inf
  integer iwork(limit)
  integer last
  integer lvl
  integer l1
  integer l2
  integer l3
  integer neval
  real result
  real work(lenw)
!
  external f
!
!  check validity of limit and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00

  if ( limit < 1 .or. lenw < limit * 4 ) go to 10
!
!  prepare call for qagie.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2

  call qagie(f,bound,inf,epsabs,epsrel,limit,result,abserr, &
    neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
  lvl = 0

10    continue

  if ( ier==6) lvl = 1
  if ( ier/=0) call xerror( 'abnormal return from  qagi', 26,ier,lvl)

  return
end
subroutine qagie(f,bound,inf,epsabs,epsrel,limit,result,abserr, &
  neval,ier,alist,blist,rlist,elist,iord,last)
!
!*******************************************************************************
!
!! QAGIE ??
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
  integer, parameter :: limexp = 50
  integer limit
!
  real abseps
  real abserr
  real alist(limit)
  real area
  real area1
  real area12
  real area2
  real a1
  real a2
  real blist(limit)
  real boun
  real bound
  real b1
  real b2
  real correc
  real defabs
  real defab1
  real defab2
  real dres
  real elist(limit)
  real epmach
  real epsabs
  real epsrel
  real erlarg
  real erlast
  real errbnd
  real errmax
  real error1
  real error2
  real erro12
  real errsum
  real ertest
  logical extrap
  real f
  integer id
  integer ier
  integer ierr
  integer ierro
  integer inf
  integer iord(limit)
  integer iroff1
  integer iroff2
  integer iroff3
  integer jupbnd
  integer k
  integer ksgn
  integer ktmin
  integer last
  logical lerr
  integer maxerr
  integer neval
  logical newflg
  logical noext
  integer nrmax
  real resabs
  real reseps
  real result
  real rlist(limit)
  real rlist2(limexp+7)
  real r1mach
  real small
  real uflow
!
!           limexp is the size of the epsilon table that can be
!           generated in ea
!
  external f
!
   epmach = epsilon ( epmach )
!
!           test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = 0.0E+00
  blist(1) = 1.0E+00
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00
  iord(1) = 0
  newflg = .true.

  if ( epsabs<=0.0E+00 .and. epsrel< max (0.5e+02*epmach,0.5e-14) ) then
    ier = 6
    go to 999
  end if
!
!  first approximation to the integral
!
!  determine the interval to be mapped onto (0,1).
!  if inf = 2 the integral is computed as i = i1+i2, where
!  i1 = integral of f over (-infinity,0),
!  i2 = integral of f over (0,+infinity).
!
  boun = bound
  if ( inf==2) boun = 0.0E+00
  call qk15i(f,boun,inf,0.0,1.0,result,abserr, defabs,resabs)
!
!  test on accuracy
!
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  dres = abs(result)
  errbnd = max (epsabs,epsrel*dres)
  if ( abserr<=1.0e+02*epmach*defabs.and.abserr>errbnd) ier = 2
  if ( limit==1) ier = 1

  if ( ier/=0.or.(abserr<=errbnd.and.abserr/=resabs).or.abserr==0.0) then
    go to 130
  end if
!
!  Initialization.
!
  uflow = 2.0E+00 * tiny ( uflow )
  lerr = .false.
  call ea(newflg,result,limexp,reseps,abseps,rlist2,ierr)
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  ktmin = 0
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ksgn = -1
  if ( dres>=(1.0-0.5e+02*epmach)*defabs) ksgn = 1
!
!  main do-loop
!
  do 90 last = 2,limit
!
!  bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1)
    call qk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2)
!
!  Improve previous approximations to integral
!  and error and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if ( defab1==error1.or.defab2==error2)go to 15

    if ( abs(rlist(maxerr)-area12)>0.1e-04*abs(area12) .or. &
      erro12<0.99*errmax) go to 10

    if ( extrap) iroff2 = iroff2+1
    if ( .not.extrap) iroff1 = iroff1+1
   10   if ( last>10.and.erro12>errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max (epsabs,epsrel*abs(area))
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2>=10.or.iroff3>=20) ier = 2
    if ( iroff2>=5) ierro = 3
!
!  Set error flag in the case that the number of subintervals equals limit.
!
    if ( last==limit) ier = 1
!
!  Set error flag in the case of bad integrand behaviour
!  at some points of the integration range.
!
    if ( max (abs(a1),abs(b2))<=(1.0+0.1e+03*epmach)* &
      (abs(a2)+0.1e+04*uflow)) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2>error1) go to 20
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 30
   20   alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
!
!           call subroutine qpsrt to maintain the descending ordering
!           in the list of error estimates and select the
!           subinterval with nrmax-th largest error estimate (to be
!           bisected next).
!
   30   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    if ( errsum<=errbnd) go to 115
    if ( ier/=0) go to 100
    if ( last==2) go to 80
    if ( noext) go to 90
    erlarg = erlarg-erlast
    if ( abs(b1-a1)>small) erlarg = erlarg+erro12
    if ( extrap) go to 40
!
!           test whether the interval to be bisected next is the
!           smallest interval.
!
    if ( abs(blist(maxerr)-alist(maxerr))>small) go to 90
    extrap = .true.
    nrmax = 2
   40   if ( ierro==3.or.erlarg<=ertest) go to 60
!
!           the smallest interval has the largest error.
!           before bisecting decrease the sum of the errors
!           over the larger intervals (erlarg) and perform
!           extrapolation.
!
    id = nrmax
    jupbnd = last
    if ( last>(2+limit/2)) then
      jupbnd = limit+3-last
    end if

    do k = id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if ( abs(blist(maxerr)-alist(maxerr))>small) go to 90
      nrmax = nrmax+1
    end do
!
!           perform extrapolation.
!
   60   call ea(newflg,area,limexp,reseps,abseps,rlist2,ierr)
    ktmin = ktmin+1

    if ( (ktmin>5).and.(abserr<0.1e-02*errsum).and.(lerr)) then
      ier = 5
    end if

    if ( (abseps>=abserr).and.(lerr)) go to 70
    ktmin = 0
    abserr = abseps
    lerr = .true.
    result = reseps
    correc = erlarg
    ertest = max (epsabs,epsrel*abs(reseps))
    if ( (abserr<=ertest).and.(lerr)) go to 100
!
!            prepare bisection of the smallest interval.
!
   70   if ( rlist2(limexp+3)==1) noext = .true.
    if ( ier==5) go to 100
    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small*0.5
    erlarg = errsum
    go to 90
   80   small = 0.375
    erlarg = errsum
    ertest = errbnd
    call ea(newflg,area,limexp,reseps,abseps,rlist2,ierr)
   90 continue
!
!           set final result and error estimate.
!
  100 if ( .not.lerr) go to 115
  if ( (ier+ierro)==0) go to 110
  if ( ierro==3) abserr = abserr+correc
  if ( ier==0) ier = 3
  if ( result/=0.0.and.area/=0.0)go to 105
  if ( abserr>errsum)go to 115
  if ( area==0.0) go to 130
  go to 110
  105 if ( abserr/abs(result)>errsum/abs(area))go to 115
!
!           test on divergence
!
  110 continue

    if ( ksgn==(-1).and. max (abs(result),abs(area))<=defabs*0.1e-01) then
      go to 130
    end if

  if ( 0.1e-01>(result/area).or.(result/area)>0.1e+03.or.errsum>abs(area)) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
  115 continue

  result = 0.0E+00
  do k = 1,last
    result = result+rlist(k)
  end do

  abserr = errsum
  130 neval = 30*last-15
  if ( inf==2) neval = 2*neval
  if ( ier>2) ier=ier-1
  999 return
end
subroutine qform(m,n,q,ldq,wa)
!
!*******************************************************************************
!
!! QFORM produces the explicit QR factorization of a matrix from its implicit form.
!
!
!  Discussion:
!
!    The QR factorization of a matrix is usually accumulated in implicit
!    form, that is, as a series of orthogonal transformations of the
!    original matrix.  This routine carries out those transformations,
!    to explicitly exhibit the factorization construced by QRFAC.
!
!  Parameters:
!
!       m is a positive integer input variable set to the number
!         of rows of a and the order of q.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       q is an m by m array. on input the full lower trapezoid in
!         the first min(m,n) columns of q contains the factored form.
!         on output q has been accumulated into a square matrix.
!
!       ldq is a positive integer input variable not less than m
!         which specifies the leading dimension of the array q.
!
!       wa is a work array of length m.
!
  integer ldq
  integer m
  integer n
!
  integer i
  integer j
  integer jm1
  integer k
  integer l
  integer minmn
  integer np1
  real q(ldq,m)
  real sum
  real temp
  real wa(m)
!
  minmn = min(m,n)

  do j = 2, minmn
    q(1:j-1,j) = 0.0E+00
  end do
!
!  Initialize remaining columns to those of the identity matrix.
!
  np1 = n + 1

  do j = n+1, m
    do i = 1, m
      q(i,j) = 0.0E+00
    end do
    q(j,j) = 1.0E+00
  end do
!
!     accumulate q from its factored form.
!
  do l = 1, minmn

    k = minmn - l + 1

    wa(k:m) = q(k:m,k)

    q(k:m,k) = 0.0E+00
    q(k,k) = 1.0E+00

    if (wa(k) /= 0.0) then

      do j = k, m
        sum = 0.0E+00
        do i = k, m
          sum = sum + q(i,j)*wa(i)
        end do
        temp = sum/wa(k)
        do i = k, m
          q(i,j) = q(i,j) - temp*wa(i)
        end do
      end do

    end if

  end do

  return
end
subroutine qk15(f,a,b,result,abserr,resabs,resasc)
!
!*******************************************************************************
!
!! QK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!           on entry
!
!              f      - real
!                     function subprogram defining the integrand
!                     function f(x). the actual name for f needs to be
!                       declared e x t e r n a l in the calling program.
!
!              a      - real: lower limit of integration
!
!              b      - real: upper limit of integration
!
!            on return
!              result - real: approximation to the integral i
!                       result is computed by applying the 15-point
!                       kronrod rule (resk) obtained by optimal addition
!                       of abscissae to the7-point gauss rule(resg).
!
!              abserr - real: estimate of the modulus of the absolute error,
!                       which should not exceed abs(i-result)
!
!              resabs - real: approximation to the integral j
!
!              resasc - real: approximation to the integral of abs(f-i/(b-a))
!                       over (a,b)
!
!  references
!
!    piessens r. et. al.,
!    "quadpack: a subroutine package for automatic integration"
!    springer, berlin 1983.
!
  real a
  real absc
  real abserr
  real b
  real centr
  real dhlgth
  real epmach
  real f
  real fc
  real fsum
  real fval1
  real fval2
  real fv1(7)
  real fv2(7)
  real hlgth
  integer j
  integer jtw
  integer jtwm1
  real resabs
  real resasc
  real resg
  real resk
  real reskh
  real result
  real r1mach
  real uflow
  real wg(4)
  real wgk(8)
  real xgk(8)
!
  external f
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 15-point kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 7-point
!                    gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 7-point gauss rule
!
!           wgk    - weights of the 15-point kronrod rule
!
!           wg     - weights of the 7-point gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
          0.9914553711208126,   0.9491079123427585, &
          0.8648644233597691,   0.7415311855993944, &
          0.5860872354676911,   0.4058451513773972, &
          0.2077849550078985,   0.0E+00              /

  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
     0.2293532201052922e-01,   0.6309209262997855e-01, &
     0.1047900103222502,   0.1406532597155259, &
     0.1690047266392679,   0.1903505780647854, &
     0.2044329400752989,   0.2094821410847278/

  data wg(1),wg(2),wg(3),wg(4)/ &
    0.1294849661688697,   0.2797053914892767, &
    0.3818300505051189,   0.4179591836734694/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 7-point gauss formula
!           resk   - result of the 15-point kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5*(a+b)
  hlgth = 0.5*(b-a)
  dhlgth = abs(hlgth)
!
!  compute the 15-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  resabs = abs(resk)

  do j=1,3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
  end do

  do j = 1,4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 0.5
  resasc = wgk(8)*abs(fc-reskh)

  do j=1,7
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do
  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00 .and. abserr /= 0.0) then
    abserr = resasc* min (1.0,(0.2e+03*abserr/resasc)**1.5)
  end if

  if ( resabs>uflow/(0.5e+02*epmach)) then
    abserr = max ((epmach*0.5e+02)*resabs,abserr)
  end if

  return
end
subroutine qk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)
!
!*******************************************************************************
!
!! QK15I ??
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
  real a
  real absc
  real absc1
  real absc2
  real abserr
  real b
  real boun
  real centr
  real dinf
  real r1mach
  real epmach
  real f
  real fc
  real fsum
  real fval1
  real fval2
  real fv1(7)
  real fv2(7)
  real hlgth
  integer inf
  integer j
  real resabs
  real resasc
  real resg
  real resk
  real reskh
  real result
  real tabsc1
  real tabsc2
  real uflow
  real wg(8)
  real wgk(8)
  real xgk(8)
!
  external f
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7), xgk(8)/ &
          0.9914553711208126,     0.9491079123427585, &
          0.8648644233597691,     0.7415311855993944, &
          0.5860872354676911,     0.4058451513773972, &
          0.2077849550078985,     0.0000000000000000/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7), wgk(8)/ &
          0.2293532201052922e-01,     0.6309209262997855e-01, &
          0.1047900103222502,     0.1406532597155259, &
          0.1690047266392679,     0.1903505780647854, &
          0.2044329400752989,     0.2094821410847278/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/ &
          0.0000000000000000,     0.1294849661688697, &
          0.0000000000000000,     0.2797053914892767, &
          0.0000000000000000,     0.3818300505051189, &
          0.0000000000000000,     0.4179591836734694/
!
  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  dinf = min(1,inf)
!
  centr = 0.5 * ( a + b )
  hlgth = 0.5 * ( b - a )
  tabsc1 = boun+dinf*(1.0-centr)/centr
  fval1 = f(tabsc1)
  if ( inf==2) fval1 = fval1+f(-tabsc1)
  fc = (fval1/centr)/centr
!
!           compute the 15-point kronrod approximation to
!           the integral, and estimate the error.
!
  resg = wg(8)*fc
  resk = wgk(8)*fc
  resabs = abs(resk)

  do j=1,7
    absc = hlgth*xgk(j)
    absc1 = centr-absc
    absc2 = centr+absc
    tabsc1 = boun+dinf*(1.0-absc1)/absc1
    tabsc2 = boun+dinf*(1.0-absc2)/absc2
    fval1 = f(tabsc1)
    fval2 = f(tabsc2)
    if ( inf==2) fval1 = fval1+f(-tabsc1)
    if ( inf==2) fval2 = fval2+f(-tabsc2)
    fval1 = (fval1/absc1)/absc1
    fval2 = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(j)*fsum
    resabs = resabs+wgk(j)*(abs(fval1)+abs(fval2))
  end do

  reskh = resk * 0.5
  resasc = wgk(8)*abs(fc-reskh)
  do j=1,7
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resasc = resasc*hlgth
  resabs = resabs*hlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00 .and. abserr /= 0.0E+00 ) then
    abserr = resasc* min (1.0,(0.2e+03*abserr/resasc)**1.5)
  end if

  if ( resabs>uflow/(0.5e+02*epmach)) then
    abserr = max ((epmach*0.5e+02)*resabs,abserr)
  end if

  return
end
subroutine qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
!
!*******************************************************************************
!
!! QPSRT ??
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
  integer last
!
  real elist(last)
  real ermax
  real errmax
  real errmin
  integer i
  integer ibeg
  integer ido
  integer iord(last)
  integer isucc
  integer j
  integer jbnd
  integer jupbn
  integer k
  integer limit
  integer maxerr
  integer nrmax
!
!  check whether the list contains more than two error estimates.
!
  if ( last<=2) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
!
!  this part of the routine is only executed
!  if, due to a difficult integrand, subdivision
!  increased the error estimate. in the normal case
!  the insert procedure should start after the
!  nrmax-th largest error estimate.
!
  errmax = elist(maxerr)
  if ( nrmax==1) go to 30
  ido = nrmax-1

  do i = 1,ido
    isucc = iord(nrmax-1)
    if ( errmax<=elist(isucc)) go to 30
    iord(nrmax) = isucc
    nrmax = nrmax-1
  end do
!
!           compute the number of elements in the list to
!           be maintained in descending order. this number
!           depends on the number of subdivisions still
!           allowed.
!
   30 continue

  jupbn = last
  if ( last>(limit/2+2)) jupbn = limit+3-last
  errmin = elist(last)
!
!           insert errmax by traversing the list top-down,
!           starting comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i = ibeg, jbnd
    isucc = iord(i)
    if ( errmax>=elist(isucc)) go to 60
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!           insert errmin by traversing the list bottom-up.
!
60 continue

  iord(i-1) = maxerr
  k = jbnd
  do j=i,jbnd
    isucc = iord(k)
    if ( errmin<elist(isucc)) go to 80
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last
!
!  set maxerr and ermax.
!
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end
subroutine qraux1(nr,n,r,i)
!
!*******************************************************************************
!
!! QRAUX1 interchanges rows i,i+1 of the upper hessenberg matrix r, columns i to n
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of matrix
! r(n,n)      <--> upper hessenberg matrix
! i            --> index of row to interchange (i<n)
!
  integer n
  integer nr
!
  integer i
  integer j
  real r(nr,n)
!
  do j = i, n
    call r_swap ( r(i,j), r(i+1,j) )
  end do

  return
end
subroutine qraux2(nr,n,r,i,a,b)
!
!*******************************************************************************
!
!! QRAUX2 pre-multiplies r by the jacobi rotation j(i,i+1,a,b)
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of matrix
! r(n,n)      <--> upper hessenberg matrix
! i            --> index of row
! a            --> scalar
! b            --> scalar
!
  integer n
  integer nr
!
  real a
  real b
  real c
  real den
  integer i
  integer j
  real r(nr,1)
  real s
  real y
  real z
!
  den=sqrt(a*a + b*b)
  c=a/den
  s=b/den

  do j=i,n
    y=r(i,j)
    z=r(i+1,j)
    r(i,j)=c*y - s*z
    r(i+1,j)=s*y + c*z
  end do

  return
end
subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,sigma,acnorm,wa)
!
!*******************************************************************************
!
!! QRFAC uses householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           t
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.
!
!  Parameters:
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a contains the matrix for
!         which the qr factorization is to be computed. on output
!         the strict upper trapezoidal part of a contains the strict
!         upper trapezoidal part of r, and the lower trapezoidal
!         part of a contains a factored form of q (the non-trivial
!         elements of the u vectors described above).
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       pivot is a logical input variable. if pivot is set .true.,
!         then column pivoting is enforced. if pivot is set .false.,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is .false., ipvt is not referenced.
!
!       lipvt is a positive integer input variable. if pivot is
!             .false., then lipvt may be as small as 1. if pivot is
!             .true., then lipvt must be at least n.
!
!       sigma is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with sigma.
!
!       wa is a work array of length n. if pivot is .false., then wa
!         can coincide with sigma.
!
  integer lipvt
  integer lda
  integer n
!
  real a(lda,n)
  real acnorm(n)
  real ajnorm
  real epsmch
  integer i
  integer ipvt(lipvt)
  integer j
  integer jp1
  integer k
  integer kmax
  integer m
  integer minmn
  logical pivot
  real r1mach
  real sigma(n)
  real snrm2
  real sum
  real temp
  real wa(n)
!
  epsmch = epsilon ( epsmch )
!
!  compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = snrm2(m,a(1,j),1)
    sigma(j) = acnorm(j)
    wa(j) = sigma(j)
    if (pivot) ipvt(j) = j
  end do
!
!  reduce a to r with householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn

    if ( pivot ) then
!
!  bring the column of largest norm into the pivot position.
!
      kmax = j
      do k = j, n
        if (sigma(k) > sigma(kmax)) kmax = k
      end do

      if ( kmax /= j) then

        do i = 1, m
          temp = a(i,j)
          a(i,j) = a(i,kmax)
          a(i,kmax) = temp
        end do

        sigma(kmax) = sigma(j)
        wa(kmax) = wa(j)
        k = ipvt(j)
        ipvt(j) = ipvt(kmax)
        ipvt(kmax) = k

      end if

    end if
!
!  Compute the householder transformation to reduce the
!  j-th column of a to a multiple of the j-th unit vector.
!
     ajnorm = snrm2(m-j+1,a(j,j),1)

     if (ajnorm /= 0.0) then

       if (a(j,j) < 0.0) then
         ajnorm = -ajnorm
       end if

       a(j:m,j) = a(j:m,j)/ajnorm
       a(j,j) = a(j,j) + 1.0E+00
!
!  apply the transformation to the remaining columns and update the norms.
!
       jp1 = j + 1

       do k = j+1, n

         sum = 0.0E+00
         do i = j, m
           sum = sum + a(i,j)*a(i,k)
         end do
         temp = sum/a(j,j)

         do i = j, m
           a(i,k) = a(i,k) - temp*a(i,j)
         end do

         if ( pivot .and. sigma(k) /= 0.0E+00 ) then

           temp = a(j,k)/sigma(k)
           sigma(k) = sigma(k)*sqrt( max (0.0,1.0-temp**2))

           if ( 0.05 * (sigma(k)/wa(k))**2 <= epsmch ) then
             sigma(k) = snrm2(m-j,a(jp1,k),1)
             wa(k) = sigma(k)
           end if

         end if

       end do

    end if

    sigma(j) = -ajnorm

  end do

  return
end
subroutine qrupdt(nr,n,a,u,v)
!
!*******************************************************************************
!
!! QRUPDT updates a QR factorization.
!
!
!  Discussion:
!
!    The routine finds an orthogonal N by N matrix Q* and an upper triangular
!    N by N matrix R* such that (Q*)(R*) = R + U*V'
!
!  Parameters:
!
!    Input, integer NR, the row dimension of the matrix.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, real A(NR,N), on input, contains the original QR
!    factorization.  On output, contains the revised factorization.
!
!    Input, real U(N), V(N), vectors that describe the rank one update
!    applied to the original matrix A.
!
  integer n
  integer nr
!
  real a(nr,1)
  integer i
  integer k
  real t1
  real t2
  real u(n)
  real v(n)
!
!  determine last non-zero in u(.)
!
  k=n

   10 continue

  if ( u(k) == 0.0E+00 .and. k > 1 ) then
    k = k - 1
    go to 10
  end if
!
! (k-1) jacobi rotations transform
!     r + u(v+) --> (r*) + (u(1)*e1)(v+)
! which is upper hessenberg
!
  if ( k > 1 ) then

    do i = k-1, 1, -1

      if ( u(i) == 0.0E+00 ) then
        call qraux1(nr,n,a,i)
        u(i)=u(i+1)
      else
        call qraux2(nr,n,a,i,u(i),-u(i+1))
        u(i)=sqrt(u(i)*u(i) + u(i+1)*u(i+1))
      end if

    end do

  end if
!
! r <-- r + (u(1)*e1)(v+)
!
  a(1,1:n) = a(1,1:n) + u(1) * v(1:n)
!
!  (k-1) jacobi rotations transform upper hessenberg r
!  to upper triangular (r*)
!
    do i = 1, k-1

      if ( a(i,i) == 0.0E+00 ) then
        call qraux1(nr,n,a,i)
      else
        t1=a(i,i)
        t2=-a(i+1,i)
        call qraux2(nr,n,a,i,t1,t2)
      end if

    end do

  return
end
function r1mach ( i )
!
!*******************************************************************************
!
!! R1MACH returns single precision machine constants.
!
!
!  Assume that single precision numbers are stored with a mantissa of T digits
!  in base B, with an exponent whose value must lie between EMIN and EMAX.  Then
!  for values of I between 1 and 5, R1MACH will return the following values:
!
!    R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!    R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!    R1MACH(3) = B**(-T), the smallest relative spacing.
!    R1MACH(4) = B**(1-T), the largest relative spacing.
!    R1MACH(5) = log10(B)
!
!  To alter this function for a particular environment, the desired set of data
!  statements should be activated by removing the C from column 1.
!
!  On rare machines a STATIC statement may need to be added.  But probably more
!  systems prohibit it that require it.
!
!  For IEEE-arithmetic machines (binary standard), the first set of constants
!  below should be appropriate.
!
!  Where possible, octal or hexadecimal constants have been used to specify the
!  constants exactly which has in some cases required the use of EQUIVALENCED
!  integer arrays.  If your compiler uses half-word integers by default
!  (sometimes called INTEGER*2), you may need to change INTEGER to INTEGER*4 or
!  otherwise instruct your compiler to use full-word integers in the next 5
!  declarations.
!
  integer diver(2)
  integer i
  integer large(2)
  integer log10(2)
  real r1mach
  integer right(2)
  real rmach(5)
  integer small(2)
!
  equivalence (rmach(1),small(1))
  equivalence (rmach(2),large(1))
  equivalence (rmach(3),right(1))
  equivalence (rmach(4),diver(1))
  equivalence (rmach(5),log10(1))
!
!  IEEE arithmetic machines, such as the ATT 3B series, Motorola 68000 based
!  machines such as the SUN 3 and ATT PC 7300, and 8087 based micros such as
!  the IBM PC and ATT 6300.
!
   data small(1) /     8388608 /
   data large(1) /  2139095039 /
   data right(1) /   864026624 /
   data diver(1) /   872415232 /
   data log10(1) /  1050288283 /
!
!  ALLIANT FX/8 UNIX Fortran compiler with the -r8 command line option.  This
!  option causes all variables declared with 'REAL' to be of type 'REAL*8' or
!  DOUBLE PRECISION.  This option does not override the 'REAL*4' declarations.
!  These R1MACH numbers below and the coresponding I1MACH are simply the DOUBLE
!  PRECISION or 'REAL*8' numbers.  If you use the -r8 your whole code (and the
!  user libraries you link with, the system libraries are taken care of
!  automagicly) must be compiled with this option.
!
!      data rmach(1) / 2.22507385850721D-308 /
!      data rmach(2) / 1.79769313486231D+308 /
!      data rmach(3) / 1.1101827117665D-16 /
!      data rmach(4) / 2.2203654423533D-16 /
!      data rmach(5) / 3.01029995663981E-1 /
!
!  AMDAHL machines.
!
!      data small(1) /    1048576 /
!      data large(1) / 2147483647 /
!      data right(1) /  990904320 /
!      data diver(1) / 1007681536 /
!      data log10(1) / 1091781651 /
!
!  BURROUGHS 1700 system.
!
!      data rmach(1) / Z400800000 /
!      data rmach(2) / Z5FFFFFFFF /
!      data rmach(3) / Z4E9800000 /
!      data rmach(4) / Z4EA800000 /
!      data rmach(5) / Z500E730E8 /
!
!  BURROUGHS 5700/6700/7700 systems.
!
!      data rmach(1) / O1771000000000000 /
!      data rmach(2) / O0777777777777777 /
!      data rmach(3) / O1311000000000000 /
!      data rmach(4) / O1301000000000000 /
!      data rmach(5) / O1157163034761675 /
!
!  CDC CYBER 170/180 series using NOS
!
!      data rmach(1) / O"00014000000000000000" /
!      data rmach(2) / O"37767777777777777777" /
!      data rmach(3) / O"16404000000000000000" /
!      data rmach(4) / O"16414000000000000000" /
!      data rmach(5) / O"17164642023241175720" /
!
!  CDC CYBER 170/180 series using NOS/VE
!
!      data rmach(1) / Z"3001800000000000" /
!      data rmach(2) / Z"4FFEFFFFFFFFFFFE" /
!      data rmach(3) / Z"3FD2800000000000" /
!      data rmach(4) / Z"3FD3800000000000" /
!      data rmach(5) / Z"3FFF9A209A84FBCF" /
!
!  CDC CYBER 200 series
!
!      data rmach(1) / X'9000400000000000' /
!      data rmach(2) / X'6FFF7FFFFFFFFFFF' /
!      data rmach(3) / X'FFA3400000000000' /
!      data rmach(4) / X'FFA4400000000000' /
!      data rmach(5) / X'FFD04D104D427DE8' /
!
!  CDC 6000/7000 series using FTN4.
!
!      data rmach(1) / 00564000000000000000B /
!      data rmach(2) / 37767777777777777776B /
!      data rmach(3) / 16414000000000000000B /
!      data rmach(4) / 16424000000000000000B /
!      data rmach(5) / 17164642023241175720B /
!
!  CDC 6000/7000 series using FTN5.
!
!      data rmach(1) / O"00564000000000000000" /
!      data rmach(2) / O"37767777777777777776" /
!      data rmach(3) / O"16414000000000000000" /
!      data rmach(4) / O"16424000000000000000" /
!      data rmach(5) / O"17164642023241175720" /
!
!  CONVEX C-1.
!
!      data rmach(1) / '00800000'X /
!      data rmach(2) / '7FFFFFFF'X /
!      data rmach(3) / '34800000'X /
!      data rmach(4) / '35000000'X /
!      data rmach(5) / '3F9A209B'X /
!
!  CONVEX C-120 (native mode) without -R8 option
!
!      data rmach(1) / 2.9387360E-39 /
!      data rmach(2) / 1.7014117E+38 /
!      data rmach(3) / 5.9604645E-08 /
!      data rmach(4) / 1.1920929E-07 /
!      data rmach(5) / 3.0102999E-01 /
!
!  CONVEX C-120 (native mode) with -R8 option
!
!      data rmach(1) / 5.562684646268007D-309 /
!      data rmach(2) / 8.988465674311577D+307 /
!      data rmach(3) / 1.110223024625157D-016 /
!      data rmach(4) / 2.220446049250313D-016 /
!      data rmach(5) / 3.010299956639812D-001 /
!
!  CONVEX C-120 (IEEE mode) without -R8 option
!
!      data rmach(1) / 1.1754945E-38 /
!      data rmach(2) / 3.4028234E+38 /
!      data rmach(3) / 5.9604645E-08 /
!      data rmach(4) / 1.1920929E-07 /
!      data rmach(5) / 3.0102999E-01 /
!
!  CONVEX C-120 (IEEE mode) with -R8 option
!
!      data rmach(1) / 2.225073858507202D-308 /
!      data rmach(2) / 1.797693134862315D+308 /
!      data rmach(3) / 1.110223024625157D-016 /
!      data rmach(4) / 2.220446049250313D-016 /
!      data rmach(5) / 3.010299956639812D-001 /
!
!  CRAY 1, 2, XMP and YMP.
!
!      data rmach(1) / 200034000000000000000B /
!      data rmach(2) / 577767777777777777776B /
!      data rmach(3) / 377224000000000000000B /
!      data rmach(4) / 377234000000000000000B /
!      data rmach(5) / 377774642023241175720B /
!
!  DATA GENERAL ECLIPSE S/200.
!  Note - It may be appropriate to include the line: STATIC RMACH(5)
!
!      data small /20K,0/
!      data large /77777K,177777K/
!      data right /35420K,0/
!      data diver /36020K,0/
!      data log10 /40423K,42023K/
!
!  ELXSI 6400, assuming REAL*4 is the default real type.
!
!      data small(1) / '00800000'X /
!      data large(1) / '7F7FFFFF'X /
!      data right(1) / '33800000'X /
!      data diver(1) / '34000000'X /
!      data log10(1) / '3E9A209B'X /
!
!  HARRIS 220
!
!      data small(1),small(2) / '20000000, '00000201 /
!      data large(1),large(2) / '37777777, '00000177 /
!      data right(1),right(2) / '20000000, '00000352 /
!      data diver(1),diver(2) / '20000000, '00000353 /
!      data log10(1),log10(2) / '23210115, '00000377 /
!
!  HARRIS SLASH 6 and SLASH 7.
!
!      data small(1),small(2) / '20000000, '00000201 /
!      data large(1),large(2) / '37777777, '00000177 /
!      data right(1),right(2) / '20000000, '00000352 /
!      data diver(1),diver(2) / '20000000, '00000353 /
!      data log10(1),log10(2) / '23210115, '00000377 /
!
!  HONEYWELL DPS 8/70 and 600/6000 series.
!
!      data rmach(1) / O402400000000 /
!      data rmach(2) / O376777777777 /
!      data rmach(3) / O714400000000 /
!      data rmach(4) / O716400000000 /
!      data rmach(5) / O776464202324 /
!
!  HP 2100, 3 word double precision with FTN4
!
!      data small(1), small(2) / 40000B,       1 /
!      data large(1), large(2) / 77777B, 177776B /
!      data right(1), right(2) / 40000B,    325B /
!      data diver(1), diver(2) / 40000B,    327B /
!      data log10(1), log10(2) / 46420B,  46777B /
!
!  HP 2100, 4 word double precision with FTN4
!
!      data small(1), small(2) / 40000B,       1 /
!      data large91), large(2) / 77777B, 177776B /
!      data right(1), right(2) / 40000B,    325B /
!      data diver(1), diver(2) / 40000B,    327B /
!      data log10(1), log10(2) / 46420B,  46777B /
!
!  HP 9000
!
!      r1mach(1) = 1.17549435E-38
!      r1mach(2) = 1.70141163E+38
!      r1mach(3) = 5.960464478E-8
!      r1mach(4) = 1.119209290E-7
!      r1mach(5) = 3.01030010E-1
!
!      data small(1) / 00040000000B /
!      data large(1) / 17677777777B /
!      data right(1) / 06340000000B /
!      data diver(1) / 06400000000B /
!      data log10(1) / 07646420233B /
!
!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!  and PERKIN ELMER (INTERDATA) 3230.
!
!      data rmach(1) / Z00100000 /
!      data rmach(2) / Z7FFFFFFF /
!      data rmach(3) / Z3B100000 /
!      data rmach(4) / Z3C100000 /
!      data rmach(5) / Z41134413 /
!
!  IBM PC - Microsoft FORTRAN
!
!      data small(1) / #00800000 /
!      data large(1) / #7F7FFFFF /
!      data right(1) / #33800000 /
!      data diver(1) / #34000000 /
!      data log10(1) / #3E9A209A /
!
!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!
!      data small(1)/ Z'00800000' /
!      data large(1)/ Z'7F7FFFFF' /
!      data right(1)/ Z'33800000' /
!      data diver(1)/ Z'34000000' /
!      data log10(1)/ Z'3E9A209A' /
!
!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!  For the INTERDATA FORTRAN VII compiler replace the Z'S specifying HEX
!  constants with Y'S.
!
!      data rmach(1) / Z'00100000' /
!      data rmach(2) / Z'7EFFFFFF' /
!      data rmach(3) / Z'3B100000' /
!      data rmach(4) / Z'3C100000' /
!      data rmach(5) / Z'41134413' /
!
!  PDP-10 (KA or KI processor).
!
!      data rmach(1) / "000400000000 /
!      data rmach(2) / "377777777777 /
!      data rmach(3) / "146400000000 /
!      data rmach(4) / "147400000000 /
!      data rmach(5) / "177464202324 /
!
!  PDP-11 FORTRANS supporting 32-bit integers (integer version).
!
!      data small(1) /    8388608 /
!      data large(1) / 2147483647 /
!      data right(1) /  880803840 /
!      data diver(1) /  889192448 /
!      data log10(1) / 1067065499 /
!
!  PDP-11 FORTRANS supporting 32-bit integers (octal version).
!
!      data rmach(1) / O00040000000 /
!      data rmach(2) / O17777777777 /
!      data rmach(3) / O06440000000 /
!      data rmach(4) / O06500000000 /
!      data rmach(5) / O07746420233 /
!
!  PDP-11 FORTRANS supporting 16-bit integers (integer version).
!
!      data small(1),small(2) /   128,     0 /
!      data large(1),large(2) / 32767,    -1 /
!      data right(1),right(2) / 13440,     0 /
!      data diver(1),diver(2) / 13568,     0 /
!      data log10(1),log10(2) / 16282,  8347 /
!
!  PDP-11 FORTRANS supporting 16-bit integers (octal version).
!
!      data small(1),small(2) / O000200, O000000 /
!      data large(1),large(2) / O077777, O177777 /
!      data right(1),right(2) / O032200, O000000 /
!      data diver(1),diver(2) / O032400, O000000 /
!      data log10(1),log10(2) / O037632, O020233 /
!
!  SEQUENT BALANCE 8000.
!
!      data small(1) / $00800000 /
!      data large(1) / $7F7FFFFF /
!      data right(1) / $33800000 /
!      data diver(1) / $34000000 /
!      data log10(1) / $3E9A209B /
!
!  SUN Microsystems UNIX F77 compiler.
!
!      data rmach(1) / 1.17549435E-38 /
!      data rmach(2) / 3.40282347E+38 /
!      data rmach(3) / 5.96016605E-08 /
!      data rmach(4) / 1.19203321E-07 /
!      data rmach(5) / 3.01030010E-01 /
!
!  SUN 3 (68881 or FPA)
!
!      data small(1) / X'00800000' /
!      data large(1) / X'7F7FFFFF' /
!      data right(1) / X'33800000' /
!      data diver(1) / X'34000000' /
!      data log10(1) / X'3E9A209B' /
!
!  UNIVAC 1100 series.
!
!      data rmach(1) / O000400000000 /
!      data rmach(2) / O377777777777 /
!      data rmach(3) / O146400000000 /
!      data rmach(4) / O147400000000 /
!      data rmach(5) / O177464202324 /
!
!  VAX/ULTRIX F77 compiler.
!
!      data small(1) /       128 /
!      data large(1) /    -32769 /
!      data right(1) /     13440 /
!      data diver(1) /     13568 /
!      data log10(1) / 547045274 /
!
!  VAX-11 with FORTRAN IV-PLUS compiler.
!
!      data rmach(1) / Z00000080 /
!      data rmach(2) / ZFFFF7FFF /
!      data rmach(3) / Z00003480 /
!      data rmach(4) / Z00003500 /
!      data rmach(5) / Z209B3F9A /
!
!  VAX/VMS version 2.2.
!
!      data rmach(1) /       '80'X /
!      data rmach(2) / 'FFFF7FFF'X /
!      data rmach(3) /     '3480'X /
!      data rmach(4) /     '3500'X /
!      data rmach(5) / '209B3F9A'X /
!
!  VAX/VMS 11/780
!
!      data small(1) / Z00000080 /
!      data large(1) / ZFFFF7FFF /
!      data right(1) / Z00003480 /
!      data diver(1) / Z00003500 /
!      data log10(1) / Z209B3F9A /
!
!  Z80 microprocessor.
!
!      data small(1), small(2) /     0,    256 /
!      data large(1), large(2) /    -1,   -129 /
!      data right(1), right(2) /     0,  26880 /
!      data diver(1), diver(2) /     0,  27136 /
!      data log10(1), log10(2) /  8347,  32538 /
!
  if ( i < 1 .or. i > 5 ) then
    write ( *, * ) ' '
    write(*,*)'R1MACH - Fatal error!'
    write(*,*)'I is out of bounds=',i
    r1mach=0.0E+00
    stop
  else
    r1mach = rmach(i)
  end if

  return
end
subroutine r1mpyq ( m, n, a, lda, v, w )
!
!*******************************************************************************
!
!! R1MPYQ is given an m by n matrix a, this subroutine computes a*q where
!     q is the product of 2*(n - 1) transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     and gv(i), gw(i) are givens rotations in the (i,n) plane which
!     eliminate elements in the i-th and n-th planes, respectively.
!     q itself is not given, rather the information to recover the
!     gv, gw rotations is supplied.
!
!  Parameters:
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a must contain the matrix
!         to be postmultiplied by the orthogonal matrix q
!         described above. on output a*q has replaced a.
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       v is an input array of length n. v(i) must contain the
!         information necessary to recover the givens rotation gv(i)
!         described above.
!
!       w is an input array of length n. w(i) must contain the
!         information necessary to recover the givens rotation gw(i)
!         described above.
!
  integer lda
  integer n
!
  real a(lda,n)
  real cos2
  integer i
  integer j
  integer m
  integer nmj
  real sin2
  real temp
  real v(n)
  real w(n)
!
  do nmj = 1, n-1

    j = n - nmj

    if (abs(v(j)) > 1.0) then
      cos2 = 1.0/v(j)
      sin2 = sqrt(1.0-cos2**2)
    else
      sin2 = v(j)
      cos2 = sqrt(1.0-sin2**2)
    end if

    do i = 1, m
      temp = cos2*a(i,j) - sin2*a(i,n)
      a(i,n) = sin2*a(i,j) + cos2*a(i,n)
      a(i,j) = temp
    end do

  end do
!
!     apply the second set of givens rotations to a.
!
  do j = 1, n-1

     if (abs(w(j)) > 1.0) then
       cos2 = 1.0/w(j)
       sin2 = sqrt(1.0-cos2**2)
     else
       sin2 = w(j)
       cos2 = sqrt(1.0-sin2**2)
     end if

     do i = 1, m
        temp = cos2*a(i,j) + sin2*a(i,n)
        a(i,n) = -sin2*a(i,j) + cos2*a(i,n)
        a(i,j) = temp
     end do

  end do

  return
end
subroutine r1updt ( m, n, s, ls, u, v, w, sing )
!
!*******************************************************************************
!
!! R1UPDT updates a QR factorization after a rank 1 update of the original matrix.
!
!
!  Discussion:
!
!    The routine is given an m by n lower trapezoidal matrix s, an m-vector u,
!    and an n-vector v, the problem is to determine an
!    orthogonal matrix q such that
!
!      (s + u*v')*q
!
!    is again lower trapezoidal.
!
!    This subroutine determines q as the product of 2*(n - 1) transformations
!
!      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!    where gv(i), gw(i) are givens rotations in the (i,n) plane
!    which eliminate elements in the i-th and n-th planes,
!    respectively. q itself is not accumulated, rather the
!    information to recover the gv, gw rotations is returned.
!
!  Parameters:
!
!       m is a positive integer input variable set to the number
!         of rows of s.
!
!       n is a positive integer input variable set to the number
!         of columns of s. n must not exceed m.
!
!       s is an array of length ls. on input s must contain the lower
!         trapezoidal matrix s stored by columns. on output s contains
!         the lower trapezoidal matrix produced as described above.
!
!       ls is a positive integer input variable not less than
!         (n*(2*m-n+1))/2.
!
!       u is an input array of length m which must contain the
!         vector u.
!
!       v is an array of length n. on input v must contain the vector
!         v. on output v(i) contains the information necessary to
!         recover the givens rotation gv(i) described above.
!
!       w is an output array of length m. w(i) contains information
!         necessary to recover the givens rotation gw(i) described
!         above.
!
!       sing is a logical output variable. sing is set .true. if any
!         of the diagonal elements of the output s are zero. otherwise
!         sing is set .false.
!
  integer ls
  integer m
  integer n
!
  real cos
  real cotan
  real giant
  integer i
  integer j
  integer jj
  integer l
  integer nmj
  real, parameter :: p25 = 0.25E+00
  real, parameter :: p5 = 0.5E+00
  real r1mach
  real s(ls)
  real sin
  logical sing
  real tan
  real tau
  real temp
  real u(m)
  real v(n)
  real w(m)
!
  giant = huge ( giant )
!
!  Initialize the diagonal element pointer.
!
  jj = (n*(2*m - n + 1))/2 - (m - n)
!
!     move the nontrivial part of the last column of s into w.
!
  l = jj
  do i = n, m
    w(i) = s(l)
    l = l + 1
  end do
!
!     rotate the vector v into a multiple of the n-th unit vector
!     in such a way that a spike is introduced into w.
!
  do nmj = 1, n-1

     j = n - nmj
     jj = jj - (m - j + 1)
     w(j) = 0.0E+00

      if ( v(j) /= 0.0E+00 ) then
!
!        determine a givens rotation which eliminates the
!        j-th element of v.
!
      if ( abs(v(n)) < abs(v(j)) ) then
        cotan = v(n)/v(j)
        sin = p5/sqrt(p25+p25*cotan**2)
        cos = sin*cotan
        tau = 1.0E+00
        if (abs(cos)*giant > 1.0) tau = 1.0/cos
      else
        tan = v(j)/v(n)
        cos = p5/sqrt(p25+p25*tan**2)
        sin = cos*tan
        tau = sin
      end if
!
!  apply the transformation to v and store the information
!  necessary to recover the givens rotation.
!
      v(n) = sin*v(j) + cos*v(n)
      v(j) = tau
!
!  apply the transformation to s and extend the spike in w.
!
      l = jj

      do i = j, m
        temp = cos*s(l) - sin*w(i)
        w(i) = sin*s(l) + cos*w(i)
        s(l) = temp
        l = l + 1
      end do

    end if

  end do
!
!     add the spike from the rank 1 update to w.
!
  w(1:m) = w(1:m) + v(n) * u(1:m)
!
!     eliminate the spike.
!
  sing = .false.

  do j = 1, n-1

     if ( w(j) /= 0.0E+00 ) then
!
!        determine a givens rotation which eliminates the
!        j-th element of the spike.
!
     if ( abs ( s(jj)) < abs ( w(j) ) ) then
        cotan = s(jj)/w(j)
        sin = p5/sqrt(p25+p25*cotan**2)
        cos = sin*cotan
        tau = 1.0E+00
        if (abs(cos)*giant > 1.0) tau = 1.0/cos
     else
        tan = w(j)/s(jj)
        cos = p5/sqrt(p25+p25*tan**2)
        sin = cos*tan
        tau = sin
     end if
!
!  apply the transformation to s and reduce the spike in w.
!
     l = jj
     do i = j, m
        temp = cos*s(l) + sin*w(i)
        w(i) = -sin*s(l) + cos*w(i)
        s(l) = temp
        l = l + 1
     end do
!
!  store the information necessary to recover the givens rotation.
!
     w(j) = tau
    end if
!
!        test for zero diagonal elements in the output s.
!
     if (s(jj) == 0.0) sing = .true.
     jj = jj + (m - j + 1)

  end do
!
!     move w back into the last column of the output s.
!
  l = jj
  do i = n, m
     s(l) = w(i)
     l = l + 1
  end do

  if (s(jj) == 0.0) sing = .true.

  return
end
function r9lgmc ( x )
!
!*******************************************************************************
!
!! R9LGMC computes the log gamma correction factor so that
!            log(gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x
!            + r9lgmc(x)
!
!  Description:
!
! compute the log gamma correction factor for x >= 10.0E+00 so that
!  log (gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x + r9lgmc(x)
!
! series for algm       on the interval  0.          to  1.00000d-02
!                                        with weighted error   3.40e-16
!                                         log weighted error  15.47
!                               significant figures required  14.39
!                                    decimal places required  15.86
!
  real, save, dimension ( 6 ) :: algmcs = (/ &
     0.166638948045186, -0.0000138494817606,  0.0000000098108256, &
    -0.0000000000180912, 0.0000000000000622, -0.0000000000000003 /)
  real csevl
  integer inits
  integer, save :: nalgm = 0
  real r1mach
  real r9lgmc
  real x
  real, save :: xbig = 0.0E+00
  real, save :: xmax = 0.0E+00
!
  if ( nalgm == 0 ) then
    nalgm = inits (algmcs, 6, r1mach(3))
    xbig = 1.0/sqrt(r1mach(3))
    xmax = exp ( min (log ( huge ( xmax ) / 12.0E+00 ), &
                     -log ( 12.0E+00 * tiny ( xmax ) )) )
  end if

  if ( x < 10.0E+00 ) then
    call xerror ( 'r9lgmc  x must be ge 10', 23, 1, 2)
  end if

  if ( x < xmax ) then
    r9lgmc = 1.0/(12.0*x)
    if (x<xbig) r9lgmc = csevl (2.0*(10./x)**2-1., algmcs, nalgm)/x
  else
    r9lgmc = 0.0E+00
    call xerror ( 'r9lgmc  x so big r9lgmc underflows', 34, 2, 1)
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
subroutine radb2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! RADB2 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
  integer ido
  integer l1
!
  real cc(ido,2,l1)
  real ch(ido,l1,2)
  integer i
  integer ic
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  ch(1,1:l1,1) = cc(1,1,1:l1) + cc(ido,2,1:l1)
  ch(1,1:l1,2) = cc(1,1,1:l1) - cc(ido,2,1:l1)

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ch(i-1,k,1) = cc(i-1,1,k) + cc(ic-1,2,k)
        tr2         = cc(i-1,1,k) - cc(ic-1,2,k)
        ch(i,k,1)   = cc(i,1,k)   - cc(ic,2,k)
        ti2         = cc(i,1,k)   + cc(ic,2,k)

        ch(i-1,k,2) = wa1(i-2) * tr2 - wa1(i-1) * ti2
        ch(i,k,2)   = wa1(i-2) * ti2 + wa1(i-1) * tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(ido,1:l1,1) =    cc(ido,1,1:l1) + cc(ido,1,1:l1)
  ch(ido,1:l1,2) = -( cc(1,2,1:l1)   + cc(1,2,1:l1) )

  return
end
subroutine radb3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! RADB3 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,3,l1)
  real ch(ido,l1,3)
  real ci2
  real ci3
  real cr2
  real cr3
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer ic
  integer k
  real, parameter :: taui =  0.866025403784439E+00
  real, parameter :: taur = -0.5E+00
  real ti2
  real tr2
  real wa1(ido)
  real wa2(ido)
!
  do k = 1, l1

    tr2 = cc(ido,2,k) + cc(ido,2,k)
    cr2 = cc(1,1,k) + taur * tr2
    ch(1,k,1) = cc(1,1,k) + tr2
    ci3 = taui * ( cc(1,3,k) + cc(1,3,k) )

    ch(1,k,2) = cr2 - ci3
    ch(1,k,3) = cr2 + ci3

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      cr2 = cc(i-1,1,k) + taur * tr2
      ch(i-1,k,1) = cc(i-1,1,k) + tr2

      ti2 = cc(i,3,k) - cc(ic,2,k)
      ci2 = cc(i,1,k) + taur * ti2
      ch(i,k,1) = cc(i,1,k) + ti2

      cr3 = taui * ( cc(i-1,3,k) - cc(ic-1,2,k) )
      ci3 = taui * ( cc(i,3,k)   + cc(ic,2,k) )

      dr2 = cr2 - ci3
      dr3 = cr2 + ci3
      di2 = ci2 + cr3
      di3 = ci2 - cr3

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3

    end do
  end do

  return
end
subroutine radb4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! RADB4 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,4,l1)
  real ch(ido,l1,4)
  real ci2
  real ci3
  real ci4
  real cr2
  real cr3
  real cr4
  integer i
  integer ic
  integer k
  real, parameter :: sqrt2 = 1.414213562373095E+00
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  do k = 1, l1

    tr1 = cc(1,1,k) - cc(ido,4,k)
    tr2 = cc(1,1,k) + cc(ido,4,k)
    tr3 = cc(ido,2,k) + cc(ido,2,k)
    tr4 = cc(1,3,k) + cc(1,3,k)

    ch(1,k,1) = tr2 + tr3
    ch(1,k,2) = tr1 - tr4
    ch(1,k,3) = tr2 - tr3
    ch(1,k,4) = tr1 + tr4

  end do

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        ti1 = cc(i,1,k) + cc(ic,4,k)
        ti2 = cc(i,1,k) - cc(ic,4,k)
        ti3 = cc(i,3,k) - cc(ic,2,k)
        tr4 = cc(i,3,k) + cc(ic,2,k)

        tr1 = cc(i-1,1,k) - cc(ic-1,4,k)
        tr2 = cc(i-1,1,k) + cc(ic-1,4,k)
        ti4 = cc(i-1,3,k) - cc(ic-1,2,k)
        tr3 = cc(i-1,3,k) + cc(ic-1,2,k)

        ch(i-1,k,1) = tr2 + tr3
        cr3         = tr2 - tr3
        ch(i,k,1)   = ti2 + ti3
        ci3         = ti2 - ti3

        cr2 = tr1 - tr4
        cr4 = tr1 + tr4
        ci2 = ti1 + ti4
        ci4 = ti1 - ti4

        ch(i-1,k,2) = wa1(i-2) * cr2 - wa1(i-1) * ci2
        ch(i,k,2)   = wa1(i-2) * ci2 + wa1(i-1) * cr2
        ch(i-1,k,3) = wa2(i-2) * cr3 - wa2(i-1) * ci3
        ch(i,k,3)   = wa2(i-2) * ci3 + wa2(i-1) * cr3
        ch(i-1,k,4) = wa3(i-2) * cr4 - wa3(i-1) * ci4
        ch(i,k,4)   = wa3(i-2) * ci4 + wa3(i-1) * cr4

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = cc(1,2,k)   + cc(1,4,k)
    ti2 = cc(1,4,k)   - cc(1,2,k)
    tr1 = cc(ido,1,k) - cc(ido,3,k)
    tr2 = cc(ido,1,k) + cc(ido,3,k)

    ch(ido,k,1) = tr2 + tr2
    ch(ido,k,2) = sqrt2 * ( tr1 - ti1 )
    ch(ido,k,3) = ti2 + ti2
    ch(ido,k,4) = -sqrt2 * ( tr1 + ti1 )

  end do

  return
end
subroutine radb5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! RADB5 is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,5,l1)
  real ch(ido,l1,5)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer ic
  integer k
  real, parameter :: ti11 =  0.951056516295154E+00
  real, parameter :: ti12 =  0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real, parameter :: tr11 =  0.309016994374947E+00
  real, parameter :: tr12 = -0.809016994374947E+00
  real tr2
  real tr3
  real tr4
  real tr5
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  do k = 1, l1

    ti5 = cc(1,3,k) + cc(1,3,k)
    ti4 = cc(1,5,k) + cc(1,5,k)
    tr2 = cc(ido,2,k) + cc(ido,2,k)
    tr3 = cc(ido,4,k) + cc(ido,4,k)

    ch(1,k,1) = cc(1,1,k) + tr2 + tr3
    cr2 = cc(1,1,k) + tr11 * tr2 + tr12 * tr3
    cr3 = cc(1,1,k) + tr12 * tr2 + tr11 * tr3
    ci5 = ti11 * ti5 + ti12 * ti4
    ci4 = ti12 * ti5 - ti11 * ti4

    ch(1,k,2) = cr2 - ci5
    ch(1,k,3) = cr3 - ci4
    ch(1,k,4) = cr3 + ci4
    ch(1,k,5) = cr2 + ci5

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      ti5 = cc(i,3,k) + cc(ic,2,k)
      ti2 = cc(i,3,k) - cc(ic,2,k)
      ti4 = cc(i,5,k) + cc(ic,4,k)
      ti3 = cc(i,5,k) - cc(ic,4,k)
      tr5 = cc(i-1,3,k) - cc(ic-1,2,k)
      tr2 = cc(i-1,3,k) + cc(ic-1,2,k)
      tr4 = cc(i-1,5,k) - cc(ic-1,4,k)
      tr3 = cc(i-1,5,k) + cc(ic-1,4,k)

      ch(i-1,k,1) = cc(i-1,1,k) + tr2 + tr3
      ch(i,k,1)   = cc(i,1,k) + ti2 + ti3

      cr2 = cc(i-1,1,k) + tr11 * tr2 + tr12 * tr3
      ci2 = cc(i,1,k)   + tr11 * ti2 + tr12 * ti3
      cr3 = cc(i-1,1,k) + tr12 * tr2 + tr11 * tr3
      ci3 = cc(i,1,k)   + tr12 * ti2 + tr11 * ti3

      cr5 = ti11 * tr5 + ti12 * tr4
      ci5 = ti11 * ti5 + ti12 * ti4
      cr4 = ti12 * tr5 - ti11 * tr4
      ci4 = ti12 * ti5 - ti11 * ti4

      dr3 = cr3 - ci4
      dr4 = cr3 + ci4
      di3 = ci3 + cr4
      di4 = ci3 - cr4
      dr5 = cr2 + ci5
      dr2 = cr2 - ci5
      di5 = ci2 - cr5
      di2 = ci2 + cr5

      ch(i-1,k,2) = wa1(i-2) * dr2 - wa1(i-1) * di2
      ch(i,k,2)   = wa1(i-2) * di2 + wa1(i-1) * dr2
      ch(i-1,k,3) = wa2(i-2) * dr3 - wa2(i-1) * di3
      ch(i,k,3)   = wa2(i-2) * di3 + wa2(i-1) * dr3
      ch(i-1,k,4) = wa3(i-2) * dr4 - wa3(i-1) * di4
      ch(i,k,4)   = wa3(i-2) * di4 + wa3(i-1) * dr4
      ch(i-1,k,5) = wa4(i-2) * dr5 - wa4(i-1) * di5
      ch(i,k,5)   = wa4(i-2) * di5 + wa4(i-1) * dr5

    end do
  end do

  return
end
subroutine radbg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! RADBG is a lower level routine used by RFFTB1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real ai1
  real ai2
  real ar1
  real ar1h
  real ar2
  real ar2h
  real arg
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  real dc2
  real dcp
  real ds2
  real dsp
  integer i
  integer ic
  integer idij
  integer ik
  integer ipph
  integer is
  integer j
  integer j2
  integer jc
  integer k
  integer l
  integer lc
  integer nbd
  real pimach
  real wa(*)
!
  arg = 2.0E+00 * pimach() / real ( ip )
  dcp = cos ( arg )
  dsp = sin ( arg )
  nbd = ( ido - 1 ) / 2
  ipph = ( ip + 1 ) / 2
  ch(1:ido,1:l1,1) = cc(1:ido,1,1:l1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    ch(1,1:l1,j) =  cc(ido,j2-2,1:l1) + cc(ido,j2-2,1:l1)
    ch(1,1:l1,jc) = cc(1,j2-1,1:l1)   + cc(1,j2-1,1:l1)
  end do

  if ( ido /= 1 ) then

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ic = ido + 2 - i
            ch(i-1,k,j)  = cc(i-1,2*j-1,k) + cc(ic-1,2*j-2,k)
            ch(i-1,k,jc) = cc(i-1,2*j-1,k) - cc(ic-1,2*j-2,k)
            ch(i,k,j)    = cc(i,2*j-1,k)   - cc(ic,2*j-2,k)
            ch(i,k,jc)   = cc(i,2*j-1,k)   + cc(ic,2*j-2,k)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ic = ido + 2 - i
          ch(i-1,1:l1,j)  = cc(i-1,2*j-1,1:l1) + cc(ic-1,2*j-2,1:l1)
          ch(i-1,1:l1,jc) = cc(i-1,2*j-1,1:l1) - cc(ic-1,2*j-2,1:l1)
          ch(i,1:l1,j)    = cc(i,2*j-1,1:l1)   - cc(ic,2*j-2,1:l1)
          ch(i,1:l1,jc)   = cc(i,2*j-1,1:l1)   + cc(ic,2*j-2,1:l1)
        end do
      end do

    end if

  end if

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      c2(ik,l)  = ch2(ik,1) + ar1 * ch2(ik,2)
      c2(ik,lc) =             ai1 * ch2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2  = dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        c2(ik,l)  = c2(ik,l)  + ar2 * ch2(ik,j)
        c2(ik,lc) = c2(ik,lc) + ai2 * ch2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + ch2(1:idl1,j)
  end do

  do j = 2, ipph
    jc = ip + 2 - j
    ch(1,1:l1,j)  = c1(1,1:l1,j) - c1(1,1:l1,jc)
    ch(1,1:l1,jc) = c1(1,1:l1,j) + c1(1,1:l1,jc)
  end do

  if ( ido /= 1 ) then

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            ch(i-1,k,j)  = c1(i-1,k,j) - c1(i,k,jc)
            ch(i-1,k,jc) = c1(i-1,k,j) + c1(i,k,jc)
            ch(i,k,j)    = c1(i,k,j)   + c1(i-1,k,jc)
            ch(i,k,jc)   = c1(i,k,j)   - c1(i-1,k,jc)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          ch(i-1,1:l1,j)  = c1(i-1,1:l1,j) - c1(i,1:l1,jc)
          ch(i-1,1:l1,jc) = c1(i-1,1:l1,j) + c1(i,1:l1,jc)
          ch(i,1:l1,j)    = c1(i,1:l1,j)   + c1(i-1,1:l1,jc)
          ch(i,1:l1,jc)   = c1(i,1:l1,j)   - c1(i-1,1:l1,jc)
        end do
      end do

    end if

  end if

  if ( ido == 1 ) then
    return
  end if

  c2(1:idl1,1) = ch2(1:idl1,1)
  c1(1,1:l1,2:ip) = ch(1,1:l1,2:ip)

  if ( nbd <= l1 ) then

    is = -ido

    do j = 2, ip
      is = is + ido
      idij = is
      do i = 3, ido, 2
        idij = idij + 2
        c1(i-1,1:l1,j) = wa(idij-1) * ch(i-1,1:l1,j) - wa(idij) * ch(i,1:l1,j)
        c1(i,1:l1,j)   = wa(idij-1) * ch(i,1:l1,j)   + wa(idij) * ch(i-1,1:l1,j)
      end do
    end do

  else

    is = -ido
    do j = 2, ip
      is = is + ido
      do k = 1, l1
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          c1(i-1,k,j) = wa(idij-1) * ch(i-1,k,j) - wa(idij) * ch(i,k,j)
          c1(i,k,j)   = wa(idij-1) * ch(i,k,j)   + wa(idij) * ch(i-1,k,j)
        end do
      end do
    end do

  end if

  return
end
subroutine radf2 ( ido, l1, cc, ch, wa1 )
!
!*******************************************************************************
!
!! RADF2 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,l1,2)
  real ch(ido,2,l1)
  integer i
  integer ic
  integer k
  real ti2
  real tr2
  real wa1(ido)
!
  ch(1,1,1:l1)   = cc(1,1:l1,1) + cc(1,1:l1,2)
  ch(ido,2,1:l1) = cc(1,1:l1,1) - cc(1,1:l1,2)

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        tr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ti2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)

        ch(i,1,k) = cc(i,k,1) + ti2
        ch(ic,2,k) = ti2 - cc(i,k,1)
        ch(i-1,1,k) = cc(i-1,k,1) + tr2
        ch(ic-1,2,k) = cc(i-1,k,1) - tr2

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  ch(1,2,1:l1) = -cc(ido,1:l1,2)
  ch(ido,1,1:l1) = cc(ido,1:l1,1)

  return
end
subroutine radf3 ( ido, l1, cc, ch, wa1, wa2 )
!
!*******************************************************************************
!
!! RADF3 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,l1,3)
  real ch(ido,3,l1)
  real ci2
  real cr2
  real di2
  real di3
  real dr2
  real dr3
  integer i
  integer ic
  integer k
  real, parameter :: taui = 0.866025403784439E+00
  real, parameter :: taur = -0.5E+00
  real ti2
  real ti3
  real tr2
  real tr3
  real wa1(ido)
  real wa2(ido)
!
  do k = 1, l1
    cr2 = cc(1,k,2) + cc(1,k,3)
    ch(1,1,k) = cc(1,k,1) + cr2
    ch(1,3,k) = taui * ( cc(1,k,3) - cc(1,k,2) )
    ch(ido,2,k) = cc(1,k,1) + taur * cr2
  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)

      cr2 = dr2 + dr3
      ci2 = di2 + di3

      ch(i-1,1,k) = cc(i-1,k,1) + cr2
      ch(i,1,k)   = cc(i,k,1) + ci2

      tr2 = cc(i-1,k,1) + taur * cr2
      ti2 = cc(i,k,1) + taur * ci2
      tr3 = taui * ( di2 - di3 )
      ti3 = taui * ( dr3 - dr2 )

      ch(i-1,3,k) = tr2 + tr3
      ch(ic-1,2,k) = tr2 - tr3
      ch(i,3,k) = ti2 + ti3
      ch(ic,2,k) = ti3 - ti2

    end do
  end do

  return
end
subroutine radf4 ( ido, l1, cc, ch, wa1, wa2, wa3 )
!
!*******************************************************************************
!
!! RADF4 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,l1,4)
  real ch(ido,4,l1)
  real ci2
  real ci3
  real ci4
  real cr2
  real cr3
  real cr4
  real, parameter :: hsqt2 = 0.7071067811865475E+00
  integer i
  integer ic
  integer k
  real ti1
  real ti2
  real ti3
  real ti4
  real tr1
  real tr2
  real tr3
  real tr4
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
!
  do k = 1, l1
    tr1 = cc(1,k,2) + cc(1,k,4)
    tr2 = cc(1,k,1) + cc(1,k,3)
    ch(1,1,k) = tr1 + tr2
    ch(ido,4,k) = tr2 - tr1
    ch(ido,2,k) = cc(1,k,1) - cc(1,k,3)
    ch(1,3,k) = cc(1,k,4) - cc(1,k,2)
  end do

  if ( ido < 2 ) then
    return
  end if

  if ( ido > 2 ) then

    do k = 1, l1
      do i = 3, ido, 2

        ic = ido + 2 - i

        cr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
        ci2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
        cr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
        ci3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
        cr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
        ci4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)

        tr1 = cr2+cr4
        tr4 = cr4-cr2
        ti1 = ci2+ci4
        ti4 = ci2-ci4
        ti2 = cc(i,k,1) + ci3
        ti3 = cc(i,k,1) - ci3
        tr2 = cc(i-1,k,1) + cr3
        tr3 = cc(i-1,k,1) - cr3

        ch(i-1,1,k)  = tr1 + tr2
        ch(ic-1,4,k) = tr2 - tr1
        ch(i,1,k)    = ti1 + ti2
        ch(ic,4,k)   = ti1 - ti2
        ch(i-1,3,k)  = ti4 + tr3
        ch(ic-1,2,k) = tr3 - ti4
        ch(i,3,k)    = tr4 + ti3
        ch(ic,2,k)   = tr4 - ti3

      end do
    end do

    if ( mod ( ido, 2 ) == 1 ) then
      return
    end if

  end if

  do k = 1, l1

    ti1 = -hsqt2 * ( cc(ido,k,2) + cc(ido,k,4) )
    tr1 =  hsqt2 * ( cc(ido,k,2) - cc(ido,k,4) )

    ch(ido,1,k) = tr1 + cc(ido,k,1)
    ch(ido,3,k) = cc(ido,k,1) - tr1

    ch(1,2,k) = ti1 - cc(ido,k,3)
    ch(1,4,k) = ti1 + cc(ido,k,3)

  end do

  return
end
subroutine radf5 ( ido, l1, cc, ch, wa1, wa2, wa3, wa4 )
!
!*******************************************************************************
!
!! RADF5 is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer ido
  integer l1
!
  real cc(ido,l1,5)
  real ch(ido,5,l1)
  real ci2
  real ci3
  real ci4
  real ci5
  real cr2
  real cr3
  real cr4
  real cr5
  real di2
  real di3
  real di4
  real di5
  real dr2
  real dr3
  real dr4
  real dr5
  integer i
  integer ic
  integer k
  real, parameter :: ti11 =  0.951056516295154E+00
  real, parameter :: ti12 =  0.587785252292473E+00
  real ti2
  real ti3
  real ti4
  real ti5
  real, parameter :: tr11 =  0.309016994374947E+00
  real, parameter :: tr12 = -0.809016994374947E+00
  real tr2
  real tr3
  real tr4
  real tr5
  real wa1(ido)
  real wa2(ido)
  real wa3(ido)
  real wa4(ido)
!
  do k = 1, l1

    cr2 = cc(1,k,5) + cc(1,k,2)
    ci5 = cc(1,k,5) - cc(1,k,2)
    cr3 = cc(1,k,4) + cc(1,k,3)
    ci4 = cc(1,k,4) - cc(1,k,3)

    ch(1,1,k)   = cc(1,k,1) + cr2 + cr3
    ch(ido,2,k) = cc(1,k,1) + tr11 * cr2 + tr12 * cr3
    ch(1,3,k)   = ti11 * ci5 + ti12 * ci4
    ch(ido,4,k) = cc(1,k,1) + tr12 * cr2 + tr11 * cr3
    ch(1,5,k)   = ti12 * ci5 - ti11 * ci4

  end do

  if ( ido == 1 ) then
    return
  end if

  do k = 1, l1
    do i = 3, ido, 2

      ic = ido + 2 - i

      dr2 = wa1(i-2) * cc(i-1,k,2) + wa1(i-1) * cc(i,k,2)
      di2 = wa1(i-2) * cc(i,k,2)   - wa1(i-1) * cc(i-1,k,2)
      dr3 = wa2(i-2) * cc(i-1,k,3) + wa2(i-1) * cc(i,k,3)
      di3 = wa2(i-2) * cc(i,k,3)   - wa2(i-1) * cc(i-1,k,3)
      dr4 = wa3(i-2) * cc(i-1,k,4) + wa3(i-1) * cc(i,k,4)
      di4 = wa3(i-2) * cc(i,k,4)   - wa3(i-1) * cc(i-1,k,4)
      dr5 = wa4(i-2) * cc(i-1,k,5) + wa4(i-1) * cc(i,k,5)
      di5 = wa4(i-2) * cc(i,k,5)   - wa4(i-1) * cc(i-1,k,5)

      cr2 = dr2 + dr5
      ci5 = dr5 - dr2
      cr5 = di2 - di5
      ci2 = di2 + di5
      cr3 = dr3 + dr4
      ci4 = dr4 - dr3
      cr4 = di3 - di4
      ci3 = di3 + di4

      ch(i-1,1,k) = cc(i-1,k,1) + cr2 + cr3
      ch(i,1,k)   = cc(i,k,1)   + ci2 + ci3

      tr2 = cc(i-1,k,1) + tr11 * cr2 + tr12 * cr3
      ti2 = cc(i,k,1)   + tr11 * ci2 + tr12 * ci3
      tr3 = cc(i-1,k,1) + tr12 * cr2 + tr11 * cr3
      ti3 = cc(i,k,1)   + tr12 * ci2 + tr11 * ci3

      tr5 = ti11 * cr5 + ti12 * cr4
      ti5 = ti11 * ci5 + ti12 * ci4
      tr4 = ti12 * cr5 - ti11 * cr4
      ti4 = ti12 * ci5 - ti11 * ci4

      ch(i-1,3,k)  = tr2 + tr5
      ch(ic-1,2,k) = tr2 - tr5
      ch(i,3,k)    = ti2 + ti5
      ch(ic,2,k)   = ti5 - ti2
      ch(i-1,5,k)  = tr3 + tr4
      ch(ic-1,4,k) = tr3 - tr4
      ch(i,5,k)    = ti3 + ti4
      ch(ic,4,k)   = ti4 - ti3

    end do
  end do

  return
end
subroutine radfg ( ido, ip, l1, idl1, cc, c1, c2, ch, ch2, wa )
!
!*******************************************************************************
!
!! RADFG is a lower level routine used by RFFTF1.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
  integer idl1
  integer ido
  integer ip
  integer l1
!
  real ai1
  real ai2
  real ar1
  real ar1h
  real ar2
  real ar2h
  real arg
  real c1(ido,l1,ip)
  real c2(idl1,ip)
  real cc(ido,ip,l1)
  real ch(ido,l1,ip)
  real ch2(idl1,ip)
  real dc2
  real dcp
  real ds2
  real dsp
  integer i
  integer ic
  integer idij
  integer ik
  integer ipph
  integer is
  integer j
  integer j2
  integer jc
  integer k
  integer l
  integer lc
  integer nbd
  real pimach
  real wa(*)
!
  arg = 2.0E+00 * pimach() / real ( ip )
  dcp = cos ( arg )
  dsp = sin ( arg )
  ipph = ( ip + 1 ) / 2
  nbd = ( ido - 1 ) / 2

  if ( ido == 1 ) then

    c2(1:idl1,1) = ch2(1:idl1,1)

  else

    ch2(1:idl1,1) = c2(1:idl1,1)
    ch(1,1:l1,2:ip) = c1(1,1:l1,2:ip)

    if ( nbd <= l1 ) then

      is = -ido
      do j = 2, ip
        is = is + ido
        idij = is
        do i = 3, ido, 2
          idij = idij + 2
          do k = 1, l1
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    else

      is = -ido

      do j = 2, ip
        is = is + ido
        do k = 1, l1
          idij = is
          do i = 3, ido, 2
            idij = idij + 2
            ch(i-1,k,j) = wa(idij-1) * c1(i-1,k,j) + wa(idij) * c1(i,k,j)
            ch(i,k,j)   = wa(idij-1) * c1(i,k,j)   - wa(idij) * c1(i-1,k,j)
          end do
        end do
      end do

    end if

    if ( nbd >= l1 ) then

      do j = 2, ipph
        jc = ip + 2 - j
        do k = 1, l1
          do i = 3, ido, 2
            c1(i-1,k,j)  = ch(i-1,k,j)  + ch(i-1,k,jc)
            c1(i-1,k,jc) = ch(i,k,j)    - ch(i,k,jc)
            c1(i,k,j)    = ch(i,k,j)    + ch(i,k,jc)
            c1(i,k,jc)   = ch(i-1,k,jc) - ch(i-1,k,j)
          end do
        end do
      end do

    else

      do j = 2, ipph
        jc = ip + 2 - j
        do i = 3, ido, 2
          c1(i-1,1:l1,j)  = ch(i-1,1:l1,j)  + ch(i-1,1:l1,jc)
          c1(i-1,1:l1,jc) = ch(i,1:l1,j)    - ch(i,1:l1,jc)
          c1(i,1:l1,j)    = ch(i,1:l1,j)    + ch(i,1:l1,jc)
          c1(i,1:l1,jc)   = ch(i-1,1:l1,jc) - ch(i-1,1:l1,j)
        end do
      end do

    end if

  end if

  do j = 2, ipph
    jc = ip + 2 - j
    c1(1,1:l1,j)  = ch(1,1:l1,j)  + ch(1,1:l1,jc)
    c1(1,1:l1,jc) = ch(1,1:l1,jc) - ch(1,1:l1,j)
  end do

  ar1 = 1.0E+00
  ai1 = 0.0E+00

  do l = 2, ipph

    lc = ip + 2 - l
    ar1h = dcp * ar1 - dsp * ai1
    ai1 =  dcp * ai1 + dsp * ar1
    ar1 = ar1h

    do ik = 1, idl1
      ch2(ik,l) = c2(ik,1) + ar1 * c2(ik,2)
      ch2(ik,lc) =           ai1 * c2(ik,ip)
    end do

    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1

    do j = 3, ipph

      jc = ip + 2 - j
      ar2h = dc2 * ar2 - ds2 * ai2
      ai2 =  dc2 * ai2 + ds2 * ar2
      ar2 = ar2h

      do ik = 1, idl1
        ch2(ik,l) =  ch2(ik,l)  + ar2 * c2(ik,j)
        ch2(ik,lc) = ch2(ik,lc) + ai2 * c2(ik,jc)
      end do

    end do

  end do

  do j = 2, ipph
    ch2(1:idl1,1) = ch2(1:idl1,1) + c2(1:idl1,j)
  end do

  cc(1:ido,1,1:l1) = ch(1:ido,1:l1,1)

  do j = 2, ipph
    jc = ip + 2 - j
    j2 = j + j
    cc(ido,j2-2,1:l1) = ch(1,1:l1,j)
    cc(1,j2-1,1:l1)   = ch(1,1:l1,jc)
  end do

  if ( ido == 1 ) then
    return
  end if

  if ( nbd >= l1 ) then

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do k = 1, l1
        do i = 3, ido, 2
          ic = ido + 2 - i
          cc(i-1,j2-1,k)  = ch(i-1,k,j) + ch(i-1,k,jc)
          cc(ic-1,j2-2,k) = ch(i-1,k,j) - ch(i-1,k,jc)
          cc(i,j2-1,k)    = ch(i,k,j)   + ch(i,k,jc)
          cc(ic,j2-2,k)   = ch(i,k,jc)  - ch(i,k,j)
        end do
      end do
    end do

  else

    do j = 2, ipph
      jc = ip + 2 - j
      j2 = j + j
      do i = 3, ido, 2
        ic = ido + 2 - i
        cc(i-1,j2-1,1:l1)  = ch(i-1,1:l1,j) + ch(i-1,1:l1,jc)
        cc(ic-1,j2-2,1:l1) = ch(i-1,1:l1,j) - ch(i-1,1:l1,jc)
        cc(i,j2-1,1:l1)    = ch(i,1:l1,j)   + ch(i,1:l1,jc)
        cc(ic,j2-2,1:l1)   = ch(i,1:l1,jc)  - ch(i,1:l1,j)
      end do
    end do

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
subroutine result(nr,n,x,f,g,a,p,itncnt,iflg,ipr)
!
!*******************************************************************************
!
!! RESULT prints information about the optimization process.
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> iterate x[k]
! f            --> function value at x[k]
! g(n)         --> gradient at x[k]
! a(n,n)       --> hessian at x[k]
! p(n)         --> step taken
! itncnt       --> iteration number k
! iflg         --> flag controlling info to print
! ipr          --> device to which to send output
!
  integer n
  integer nr
!
  real a(nr,n)
  real f
  real g(n)
  integer i
  integer iflg
  integer ipr
  integer itncnt
  integer j
  real p(n)
  real x(n)
!
  write ( ipr, 903 ) itncnt

  if ( iflg /= 0 ) then
    write ( ipr, * ) ' result       step'
    write ( ipr,905) p(1:n)
  end if

  write ( ipr, * ) ' result       x(k)'
  write ( ipr, 905) x(1:n)
  write ( ipr, * ) ' result     function at x(k)'
  write ( ipr, 905) f
  write ( ipr, * ) ' result       gradient at x(k)'
  write ( ipr, 905) g(1:n)

  if ( iflg /= 0 ) then

    write ( ipr, * ) ' result       hessian at x(k)'
    do i = 1, n
      write ( ipr, 900) i
      write ( ipr, 902) (a(i,j),j=1,i)
    end do

  end if

  return

  900 format(' result     row',i5)
  902 format(' result       ',5(2x,e20.13))
  903 format(/'0result    iterate k=',i5)
  905 format(' result               ',5(2x,e20.13) )
end
subroutine rfftb ( n, r, wsave )
!
!*******************************************************************************
!
!! RFFTB computes a real periodic sequence from its Fourier coefficients.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier synthesis.
!
!    The transform is unnormalized.  A call to RFFTF followed by a call to
!    RFFTB will multiply the input sequence by N.
!
!    If N is even, the transform is defined by:
!
!      R_out(I) = R_in(1) + (-1)**(I-1) * R_in(N) + sum ( 2 <= K <= N/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    If N is odd, the transform is defined by:
!
!      R_out(I) = R_in(1) + sum ( 2 <= K <= (N+1)/2 )
!
!        + 2 * R_in(2*K-2) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!        - 2 * R_in(2*K-1) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(2*N+15), a work array.  The WSAVE array must be
!    initialized by calling RFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
  integer n
!
  real r(n)
  real wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rfftb1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rfftb1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! RFFTB1 is a lower level routine used by RFFTB.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Input/output, real C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real CH(N).
!
!    Input, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  real c(n)
  real ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer l1
  integer l2
  integer na
  integer nf
  real wa(n)
!
  nf = ifac(2)
  na = 0
  l1 = 1
  iw = 1

  do k1 = 1, nf

    ip = ifac(k1+2)
    l2 = ip * l1
    ido = n / l2
    idl1 = ido * l1

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radb4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radb4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

      na = 1 - na

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radb2 ( ido, l1, c, ch, wa(iw) )
      else
        call radb2 ( ido, l1, ch, c, wa(iw) )
      end if

      na = 1 - na

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radb3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radb3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

      na = 1 - na

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radb5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radb5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

      na = 1 - na

    else

      if ( na == 0 ) then
        call radbg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
      else
        call radbg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
      end if

      if ( ido == 1 ) then
        na = 1 - na
      end if

    end if

    l1 = l2
    iw = iw + ( ip - 1 ) * ido

  end do

  if ( na /= 0 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine rfftf ( n, r, wsave )
!
!*******************************************************************************
!
!! RFFTF computes the Fourier coefficients of a real periodic sequence.
!
!
!  Discussion:
!
!    This process is sometimes called Fourier analysis.
!
!    The transform is unnormalized.  A call to RFFTF followed by a call
!    to RFFTB will multiply the input sequence by N.
!
!    The transform is defined by:
!
!      R_out(1) = sum ( 1 <= I <= N ) R_in(I)
!
!    Letting L = (N+1)/2, then for K = 2,...,L
!
!      R_out(2*K-2) = sum ( 1 <= I <= N )
!
!        R_in(I) * cos ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!      R_out(2*K-1) = sum ( 1 <= I <= N )
!
!        -R_in(I) * sin ( ( K - 1 ) * ( I - 1 ) * 2 * PI / N )
!
!    And, if N is even, then:
!
!      R_out(N) = sum ( 1 <= I <= N ) (-1)**(I-1) * R_in(I)
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real R(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(2*N+15), a work array.  The WSAVE array must be
!    initialized by calling RFFTI.  A different WSAVE array must be used
!    for each different value of N.
!
  integer n
!
  real r(n)
  real wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rfftf1 ( n, r, wsave(1), wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rfftf1 ( n, c, ch, wa, ifac )
!
!*******************************************************************************
!
!! RFFTF1 is a lower level routine used by RFFTF and SINT.
!
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Input/output, real C(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real CH(N).
!
!    Input, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  real c(n)
  real ch(n)
  integer idl1
  integer ido
  integer ifac(15)
  integer ip
  integer iw
  integer ix2
  integer ix3
  integer ix4
  integer k1
  integer kh
  integer l1
  integer l2
  integer na
  integer nf
  real wa(n)
!
  nf = ifac(2)
  na = 1
  l2 = n
  iw = n

  do k1 = 1, nf

    kh = nf - k1
    ip = ifac(kh+3)
    l1 = l2 / ip
    ido = n / l2
    idl1 = ido * l1
    iw = iw - ( ip - 1 ) * ido
    na = 1 - na

    if ( ip == 4 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido

      if ( na == 0 ) then
        call radf4 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3) )
      else
        call radf4 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3) )
      end if

    else if ( ip == 2 ) then

      if ( na == 0 ) then
        call radf2 ( ido, l1, c, ch, wa(iw) )
      else
        call radf2 ( ido, l1, ch, c, wa(iw) )
      end if

    else if ( ip == 3 ) then

      ix2 = iw + ido

      if ( na == 0 ) then
        call radf3 ( ido, l1, c, ch, wa(iw), wa(ix2) )
      else
        call radf3 ( ido, l1, ch, c, wa(iw), wa(ix2) )
      end if

    else if ( ip == 5 ) then

      ix2 = iw + ido
      ix3 = ix2 + ido
      ix4 = ix3 + ido

      if ( na == 0 ) then
        call radf5 ( ido, l1, c, ch, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      else
        call radf5 ( ido, l1, ch, c, wa(iw), wa(ix2), wa(ix3), wa(ix4) )
      end if

    else

      if ( ido == 1 ) then
        na = 1 - na
      end if

      if ( na == 0 ) then
        call radfg ( ido, ip, l1, idl1, c, c, c, ch, ch, wa(iw) )
        na = 1
      else
        call radfg ( ido, ip, l1, idl1, ch, ch, ch, c, c, wa(iw) )
        na = 0
      end if

    end if

    l2 = l1

  end do

  if ( na /= 1 ) then
    c(1:n) = ch(1:n)
  end if

  return
end
subroutine rffti ( n, wsave )
!
!*******************************************************************************
!
!! RFFTI initializes WSAVE, used in RFFTF and RFFTB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Output, real WSAVE(2*N+15), contains data, dependent on the value
!    of N, which is necessary for the RFFTF and RFFTB routines.
!
  integer n
!
  real wsave(2*n+15)
!
  if ( n <= 1 ) then
    return
  end if

  call rffti1 ( n, wsave(n+1), wsave(2*n+1) )

  return
end
subroutine rffti1 ( n, wa, ifac )
!
!*******************************************************************************
!
!! RFFTI1 is a lower level routine used by RFFTI.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input, real WA(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  real arg
  real argh
  real argld
  real fi
  integer i
  integer ib
  integer ido
  integer ifac(15)
  integer ii
  integer ip
  integer is
  integer j
  integer k1
  integer l1
  integer l2
  integer ld
  integer nf
  real pimach
  real wa(n)
!
  call i_factor ( n, ifac )

  nf = ifac(2)

  argh = 2.0E+00 * pimach() / real ( n )
  is = 0
  l1 = 1

  do k1 = 1, nf-1

    ip = ifac(k1+2)
    ld = 0
    l2 = l1 * ip
    ido = n / l2

    do j = 1, ip-1

      ld = ld + l1
      i = is
      argld = real ( ld ) * argh
      fi = 0.0E+00

      do ii = 3, ido, 2
        i = i + 2
        fi = fi + 1.0E+00
        arg = fi * argld
        wa(i-1) = cos ( arg )
        wa(i) = sin ( arg )
      end do

      is = is + ido

    end do

    l1 = l2

  end do

  return
end
function rnor ( )
!
!*******************************************************************************
!
!! RNOR generates normal random numbers.
!
!
!  Description:
!
!    RNOR generates normal random numbers with zero mean and
!    unit standard deviation, often denoted n(0,1).
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    Marsaglia and Tsang,
!    A fast, easily implemented method for sampling from decreasing or
!      symmetric unimodal density functions,
!    SIAM Journal of Scientific and Statistical Computing, 1983.
!
!   use
!       first time....
!                   z = rstart(iseed)
!                     here iseed is any  n o n - z e r o  integer.
!                     this causes initialization of the program.
!                     rstart returns a real (single precision) echo of iseed.
!
!       subsequent times...
!                   z = rnor()
!                     causes the next real (single precision) random number
!                           to be returned as z.
!
!                 typical usage
!
!                    real rstart,rnor,z
!                    integer iseed,i
!                    iseed = 305
!                    z = rstart(iseed)
!                    do i = 1,10
!                       z = rnor()
!                       write(*,*) z
!                    end do
!                  end
!
!
  real, parameter :: aa = 12.37586
  real, parameter :: b = 0.4878992
  real, parameter :: c = 12.67706
  real, save :: c1 = 0.9689279
  real, save :: c2 = 1.301198
  integer ia
  integer ib
  integer ic
  integer id
  integer, save :: ii = 17
  integer iii
  integer iseed
  integer j
  integer, save :: jj = 5
  integer jjj
  real, save :: pc = 0.01958303
  real rnor
  real rstart
  real s
  real t
  real, save, dimension ( 17 ) :: u
  real un
  real v(65)
  real vni
  real x
  real, save :: xn = 2.776994
  real y
!
  data v/ .3409450, .4573146, .5397793, .6062427, .6631691 &
     , .7136975, .7596125, .8020356, .8417227, .8792102, .9148948 &
     , .9490791, .9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917 &
     ,1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635 &
     ,1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929 &
     ,1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454 &
     ,1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422 &
     ,1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166 &
     ,1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713 &
     ,2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117 &
     ,2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943/
!
!  Load data array in case user forgets to initialize.
!  This array is the result of calling uni 100000 times with seed 305.
!
  data u/ &
     0.8668672834288,  0.3697986366357,  0.8008968294805, &
     0.4173889774680,  0.8254561579836,  0.9640965269077, &
     0.4508667414265,  0.6451309529668,  0.1645456024730, &
     0.2787901807898,  0.06761531340295, 0.9663226330820, &
     0.01963343943798, 0.02947398211399, 0.1636231515294, &
     0.3976343250467,  0.2631008574685/
!
!  fast part...
!
!  basic generator is fibonacci
!
  un = u(ii)-u(jj)
  if ( un<0.0) un = un+1.
  u(ii) = un
!
!  u(ii) and un are uniform on [0,1)
!  vni is uniform on [-1,1)
!
  vni = un + un -1.
  ii = ii-1
  if ( ii==0)ii = 17
  jj = jj-1
  if ( jj==0)jj = 17
!        int(un(ii)*128) in range [0,127],  j is in range [1,64]
  j = mod(int(u(ii)*128),64)+1
!
!        pick sign as vni is positive or negative
!
  rnor = vni*v(j+1)
  if ( abs(rnor)<=v(j))return
!
! slow part; aa is a*f(0)
!
  x = (abs(rnor)-v(j))/(v(j+1)-v(j))
!
!          y is uniform on [0,1)
!
  y = u(ii)-u(jj)
  if ( y<0.0) y = y+1.
  u(ii) = y
  ii = ii-1
  if ( ii==0)ii = 17
  jj = jj-1
  if ( jj==0)jj = 17

  s = x+y
  if ( s>c2)go to 11
  if ( s<=c1)return
  if ( y>c-aa*exp(-.5*(b-b*x)**2))go to 11
  if ( exp(-.5*v(j+1)**2)+y*pc/v(j+1)<=exp(-.5*rnor**2))return
!
! tail part; .3601016 is 1./xn
!       y is uniform on [0,1)
!
   22 y = u(ii)-u(jj)
  if ( y<=0.0) y = y+1.
  u(ii) = y
  ii = ii-1
  if ( ii==0)ii = 17
  jj = jj-1
  if ( jj==0)jj = 17

  x = 0.3601016*log(y)
!
!       y is uniform on [0,1)
!
  y = u(ii)-u(jj)
  if ( y<=0.0) y = y+1.
  u(ii) = y
  ii = ii-1
  if ( ii==0)ii = 17
  jj = jj-1
  if ( jj==0)jj = 17
  if (  -2.*log(y)<=x**2 )go to 22
  rnor = sign(xn-x,rnor)
  return
   11 rnor = sign(b-b*x,rnor)
  return
!
!  fill
!
  entry rstart(iseed)
  if ( iseed/=0) then
!
!  generate random bit pattern in array based on given seed
!
    ii = 17
    jj = 5
    ia = mod(abs(iseed),32707)
    ib = 1111
    ic = 1947
    do iii = 1,17
      s = 0.0E+00
      t = .50
!
!             do for each of the bits of mantissa of word
!             loop  over 64 bits, enough for all known machines
!                   in single precision
!
      do jjj = 1,64

        id = ic-ia
        if ( id < 0 ) then
          id = id+32707
          s = s+t
        end if

        ia = ib
        ib = ic
        ic = id
        t = .5*t

      end do

      u(iii) = s

    end do

  end if
!
!  return floating echo of iseed.
!
  rstart = iseed
  return
end
subroutine rsftb ( n, r, azero, a, b )
!
!*******************************************************************************
!
!! RSFTB computes a "slow" backward Fourier transform of real data.
!
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Output, real R(N), the reconstructed data sequence.
!
!    Input, real AZERO, the constant Fourier coefficient.
!
!    Input, real A(N/2), B(N/2), the Fourier coefficients.
!
  integer n
!
  real a(n/2)
  real azero
  real b(n/2)
  integer i
  integer k
  real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510E+00
  real r(n)
  real theta
!
  r(1:n) = azero
  do i = 1, n
    do k = 1, n/2
      theta = real ( k * ( i - 1 ) * 2 ) * pi / real ( n )
      r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
    end do
  end do

  return
end
subroutine rsftf ( n, r, azero, a, b )
!
!*******************************************************************************
!
!! RSFTF computes a "slow" forward Fourier transform of real data.
!
!
!  Modified:
!
!    13 March 2001
!
!  Parameters:
!
!    Input, integer N, the number of data values.
!
!    Input, real R(N), the data to be transformed.
!
!    Output, real AZERO, = sum ( 1 <= I <= N ) R(I) / N.
!
!    Output, real A(N/2), B(N/2), the Fourier coefficients.
!
  integer n
!
  real a(1:n/2)
  real azero
  real b(1:n/2)
  integer i
  integer j
  real, parameter :: pi = 3.14159265358979323846264338327950288419716939937510E+00
  real r(n)
  real theta
!
  azero = sum ( r(1:n) ) / real ( n )

  do i = 1, n / 2

    a(i) = 0.0E+00
    b(i) = 0.0E+00

    do j = 1, n
      theta = real ( 2 * i * ( j - 1 ) ) * pi / real ( n )
      a(i) = a(i) + r(j) * cos ( theta )
      b(i) = b(i) + r(j) * sin ( theta )
    end do

    a(i) = a(i) / real ( n )
    b(i) = b(i) / real ( n )

    if ( i /= ( n / 2 ) ) then
      a(i) = 2.0E+00 * a(i)
      b(i) = 2.0E+00 * b(i)
    end if

  end do

  return
end
function runge ( x )
!
!*******************************************************************************
!
!! RUNGE evaluates Runge's function.
!
!
!  Modified:
!
!    11 November 1999
!
!  Formula:
!
!    RUNGE(X) = 1 / ( 1 + 25 * X**2 )
!
!  Discussion:
!
!    Runge's function is a common test case for interpolation
!    and approximation, over the interval [-1,1].
!
!  Parameters:
!
!    Input, real X, the argument of Runge's function.
!
!    Output, real RUNGE, the value of Runge's function.
!
  real runge
  real x
!
  runge = 1.0E+00 / ( 1.0E+00 + 25.0E+00 * x**2 )

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
function samax ( n, x, incx )
!
!*******************************************************************************
!
!! SAMAX returns the maximum absolute value of the entries in a vector.
!
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real SAMAX, the maximum absolute value of an element of X.
!
  integer i
  integer incx
  integer ix
  integer n
  real samax
  real x(*)
!
  if ( n <= 0 ) then

    samax = 0.0E+00

  else if ( n == 1 ) then

    samax = abs ( x(1) )

  else if ( incx == 1 ) then

    samax = abs ( x(1) )

    do i = 2, n
      if ( abs ( x(i) ) > samax ) then
        samax = abs ( x(i) )
      end if
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    samax = abs ( x(ix) )
    ix = ix + incx

    do i = 2, n
      if ( abs ( x(ix) ) > samax ) then
        samax = abs ( x(ix) )
      end if
      ix = ix + incx
    end do

  end if

  return
end
function sasum ( n, x, incx )
!
!*******************************************************************************
!
!! SASUM sums the absolute values of the entries of a vector.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(*), the vector to be examined.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real SASUM, the sum of the absolute values of the vector.
!
  integer i
  integer incx
  integer ix
  integer m
  integer n
  real sasum
  real stemp
  real x(*)
!
  stemp = 0.0E+00

  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 6 )

    do i = 1, m
      stemp = stemp + abs ( x(i) )
    end do

    do i = m+1, n, 6
      stemp = stemp + abs ( x(i)   ) + abs ( x(i+1) ) + abs ( x(i+2) ) &
                    + abs ( x(i+3) ) + abs ( x(i+4) ) + abs ( x(i+5) )
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      stemp = stemp + abs ( x(ix) )
      ix = ix + incx
    end do

  end if

  sasum = stemp

  return
end
subroutine saxpy ( n, sa, x, incx, y, incy )
!
!*******************************************************************************
!
!! SAXPY adds a constant times one vector to another.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real SA, the multiplier.
!
!    Input, real X(*), the vector to be scaled and added to Y.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, real Y(*), the vector to which a multiple of X is to
!    be added.
!
!    Input, integer INCY, the increment between successive entries of Y.
!
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real sa
  real x(*)
  real y(*)
!
  if ( n <= 0 ) then

  else if ( sa == 0.0E+00 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    y(1:n) = y(1:n) + sa * x(1:n)

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      y(iy) = y(iy) + sa * x(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine sdcor (dfdy,el,fa,h,impl,ipvt,matdim,miter,ml,mu,n, &
  nde,nq,t,users,y,yh,ywt,evalfa,save1,save2,a,d,jstate)
!
!*******************************************************************************
!
!! SDCOR computes corrections to the y array of SDRIV3.
!
!
!  in the case of functional iteration, update y directly from the
!  result of the last call to f.
!  in the case of the chord method, compute the corrector error and
!  solve the linear system with that as right hand side and dfdy as
!  coefficient matrix, using the lu decomposition if miter is 1, 2, 4,
!  or 5.
!
  integer matdim
  integer n
!
  real a(matdim,*)
  real d
  real dfdy(matdim,*)
  real el(13,12)
  logical evalfa
  real h
  integer i
  integer i1
  integer i2
  integer i3
  integer iflag
  integer impl
  integer ipvt(*)
  integer j
  integer jstate
  integer miter
  integer ml
  integer mu
  integer mw
  integer nde
  integer nq
  real save1(*)
  real save2(*)
  real snrm2
  real t
  real y(*)
  real yh(n,*)
  real ywt(*)
!
  external fa
  external users
!
  if (miter == 0) then

    save1(1:n) = (h*save2(1:n) - yh(1:n,2) - save1(1:n))/ywt(1:n)

    d = snrm2(n, save1, 1)/sqrt(real(n))

    save1(1:n) = h * save2(1:n) - yh(1:n,2)

  else if (miter == 1 .or. miter == 2) then

    if (impl == 0) then

      save2(1:n) = h * save2(1:n) - yh(1:n,2) - save1(1:n)

    else if (impl == 1) then

      if (evalfa) then
        call fa (n, t, y, a, matdim, ml, mu, nde)
        if (n == 0) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h * save2(1:n)

      do j = 1,n
        save2(1:n) = save2(1:n) - a(1:n,j)*(yh(j,2) + save1(j))
      end do

    else if (impl == 2) then

      if (evalfa) then
        call fa (n, t, y, a, matdim, ml, mu, nde)
        if (n == 0) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h*save2(1:n) - a(1:n,1)*(yh(i,2) + save1(1:n))

    end if

    call sgesl (dfdy, matdim, n, ipvt, save2, 0)

    save1(1:n) = save1(1:n) + save2(1:n)
    save2(1:n) = save2(1:n)/ywt(1:n)

    d = snrm2(n, save2, 1) / sqrt(real(n))

  else if (miter == 4 .or. miter == 5) then

    if (impl == 0) then

      save2(1:n) = h * save2(1:n) - yh(1:n,2) - save1(1:n)

    else if (impl == 1) then
      if (evalfa) then
        call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
        if (n == 0) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h * save2(1:n)

      mw = ml + 1 + mu

      do j = 1,n
        i1 = max(ml+1, mw+1-j)
        i2 = min(mw+n-j, mw+ml)
        do i = i1,i2
          i3 = i + j - mw
          save2(i3) = save2(i3) - a(i,j)*(yh(j,2) + save1(j))
        end do
      end do

    else if (impl == 2) then

      if (evalfa) then
        call fa (n, t, y, a, matdim, ml, mu, nde)
        if (n == 0) then
          jstate = 9
          return
        end if
      else
        evalfa = .true.
      end if

      save2(1:n) = h * save2(1:n) - a(1:n,1)*(yh(1:n,2) + save1(1:n))

    end if

    call sgbsl (dfdy, matdim, n, ml, mu, ipvt, save2, 0)

    save1(1:n) = save1(1:n) + save2(1:n)
    save2(1:n) = save2(1:n)/ywt(1:n)

    d = snrm2(n, save2, 1)/sqrt(real(n))

  else if (miter == 3) then

    iflag = 2

    call users (y, yh(1,2), ywt, save1, save2, t, h, el(1,nq), impl, &
      n, nde, iflag)

    if (n == 0) then
      jstate = 10
      return
    end if

    save1(1:n) = save1(1:n) + save2(1:n)
    save2(1:n) = save2(1:n) / ywt(1:n)

    d = snrm2(n, save2, 1) / sqrt(real(n))

  end if
end
subroutine sdcst ( maxord, mint, iswflg, el, tq )
!
!*******************************************************************************
!
!! SDCST sets coefficients used by the core integrator SDSTP.
!
!
!  the array el determines the basic method.
!  the array tq is involved in adjusting the step size in relation
!  to truncation error.  el and tq depend upon mint, and are calculated
!  for orders 1 to maxord(<= 12).  for each order nq, the coefficients
!  el are calculated from the generating polynomial:
!    l(t) = el(1,nq) + el(2,nq)*t + ... + el(nq+1,nq)*t**nq.
!  for the implicit adams methods, l(t) is given by
!    dl/dt = (1+t)*(2+t)* ... *(nq-1+t)/k,   l(-1) = 0,
!    where      k = factorial(nq-1).
!  for the gear methods,
!    l(t) = (1+t)*(2+t)* ... *(nq+t)/k,
!    where      k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!  for each order nq, there are three components of tq.
!
  real el(13,12)
  real factrl(12)
  real gamma(14)
  integer i
  integer iswflg
  integer j
  integer maxord
  integer mint
  integer mxrd
  real sum
  real tq(3,12)
!
  factrl(1) = 1.0E+00
  do i = 2, maxord
    factrl(i) = real(i)*factrl(i-1)
  end do
!
!  compute adams coefficients
!
  if (mint == 1) then

    gamma(1) = 1.0E+00
    do i = 1,maxord+1
      sum = 0.0E+00
      do j = 1,i
        sum = sum - gamma(j)/real(i-j+2)
      end do
      gamma(i+1) = sum
    end do

    el(1,1) = 1.0E+00
    el(2,1) = 1.0E+00
    el(2,2) = 1.0E+00
    el(3,2) = 1.0E+00

    do j = 3,maxord
      el(2,j) = factrl(j-1)
      do i = 3,j
        el(i,j) = real(j-1)*el(i,j-1) + el(i-1,j-1)
      end do
      el(j+1,j) = 1.0E+00
    end do

    do j = 2, maxord
      el(1,j) = el(1,j-1) + gamma(j)
      el(2,j) = 1.0E+00
      do i = 3,j+1
        el(i,j) = el(i,j)/(real(i-1)*factrl(j-1))
      end do
    end do

    do j = 1,maxord
      tq(1,j) = -1.0E+00 / (factrl(j)*gamma(j))
      tq(2,j) = -1.0E+00 / gamma(j+1)
      tq(3,j) = -1.0E+00 / gamma(j+2)
    end do
!
!  compute gear coefficients
!
  else if (mint == 2) then

    el(1,1) = 1.0E+00
    el(2,1) = 1.0E+00
    do j = 2,maxord
      el(1,j) = factrl(j)
      do i = 2,j
        el(i,j) = real(j)*el(i,j-1) + el(i-1,j-1)
      end do
      el(j+1,j) = 1.0E+00
    end do

    sum = 1.0E+00
    do j = 2,maxord
      sum = sum + 1.0/real(j)
      do i = 1, j+1
        el(i,j) = el(i,j)/(factrl(j)*sum)
      end do
    end do

    do j = 1,maxord
      if (j > 1) tq(1,j) = 1.0/factrl(j-1)
      tq(2,j) = real(j+1)/el(1,j)
      tq(3,j) = real(j+2)/el(1,j)
    end do

  end if
!
!  compute constants used in the stiffness test.
!  these are the ratio of tq(2,nq) for the gear
!  methods to those for the adams methods.
!
  if (iswflg == 3) then
    mxrd = min(maxord, 5)
    if (mint == 2) then
      gamma(1) = 1.0E+00
      do i = 1,mxrd
        sum = 0.0E+00
        do j = 1,i
          sum = sum - gamma(j)/real(i-j+2)
        end do
        gamma(i+1) = sum
      end do
    end if

    sum = 1.0E+00
    do i = 2,mxrd
      sum = sum + 1.0/real(i)
      el(1+i,1) = -real(i+1)*sum*gamma(i+1)
    end do

  end if
end
subroutine sdntl(eps,f,fa,hmax,hold,impl,jtask,matdim,maxord, &
  mint,miter,ml,mu,n,nde,save1,t,uround,users,y,ywt,h,mntold, &
  mtrold,nfe,rc,yh,a,convrg,el,fac,ier,ipvt,nq,nwait,rh,rmax, &
  save2,tq,trend,iswflg,jstate)
!
!*******************************************************************************
!
!! SDNTL sets parameters for SDSTP.
!
!
!  SDNTL is called on the first call to sdstp, on an internal restart, or
!  when the user has altered mint, miter, and/or h.
!
!  on the first call, the order is set to 1 and the initial derivatives
!  are calculated.  rmax is the maximum ratio by which h can be
!  increased in one step.  it is initially rminit to compensate
!  for the small initial h, but then is normally equal to rmnorm.
!  if a failure occurs (in corrector convergence or error test), rmax
!  is set at rmfail for the next increase.
!  if the caller has changed mint, or if jtask = 0, sdcst is called
!  to set the coefficients of the method.  if the caller has changed h,
!  yh must be rescaled.  if h or mint has been changed, nwait is
!  reset to nq + 2 to prevent further increases in h for that many
!  steps.  also, rc is reset.  rc is the ratio of new to old values of
!  the coefficient l(0)*h.  if the caller has changed miter, rc is
!  set to 0 to force the partials to be updated, if partials are used.
!
  integer matdim
  integer n
!
  real a(matdim,*)
  logical convrg
  real el(13,12)
  real eps
  real fac(*)
  real h
  real hmax
  real hold
  integer i
  logical ier
  integer iflag
  integer impl
  integer info
  integer ipvt(*)
  integer iswflg
  integer jstate
  integer jtask
  integer maxord
  integer mint
  integer miter
  integer ml
  integer mntold
  integer mtrold
  integer mu
  integer nde
  integer nfe
  integer nq
  integer nwait
  real oldl0
  real rc
  real rh
  real rmax
  real, parameter :: rminit = 10000.0E+00
  real save1(*)
  real save2(*)
  real smax
  real smin
  real snrm2
  real sum
  real sum0
  real t
  real tq(3,12)
  real trend
  real uround
  real y(*)
  real yh(n,*)
  real ywt(*)
!
  external f
  external fa
  external users
!
  ier = .false.
  if (jtask >= 0) then
    if (jtask == 0) then
      call sdcst (maxord, mint, iswflg,  el, tq)
      rmax = rminit
    end if
    rc = 0.0E+00
    convrg = .false.
    trend = 1.0E+00
    nq = 1
    nwait = 3

    call f ( n, t, y, save2 )

    if (n == 0) then
      jstate = 6
      return
    end if

    nfe = nfe + 1

    if (impl /= 0) then
      if (miter == 3) then
        iflag = 0
        call users (y, yh, ywt, save1, save2, t, h, el, impl, n, nde, iflag)
        if (n == 0) then
          jstate = 10
          return
        end if
      else if (impl == 1) then
        if (miter == 1 .or. miter == 2) then
          call fa (n, t, y, a, matdim, ml, mu, nde)
          if (n == 0) then
            jstate = 9
            return
          end if
          call sgefa (a, matdim, n, ipvt, info)
          if (info /= 0) then
            ier = .true.
            return
          end if
          call sgesl (a, matdim, n, ipvt, save2, 0)
        else if (miter == 4 .or. miter == 5) then
          call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)
          if (n == 0) then
            jstate = 9
            return
          end if
          call sgbfa (a, matdim, n, ml, mu, ipvt, info)
          if (info /= 0) then
            ier = .true.
            return
          end if
          call sgbsl (a, matdim, n, ml, mu, ipvt, save2, 0)
        end if

      else if (impl == 2) then

        call fa (n, t, y, a, matdim, ml, mu, nde)
        if (n == 0) then
          jstate = 9
          return
        end if

        do i = 1, nde

          if (a(i,1) == 0.0) then
            ier = .true.
            return
          else
            save2(i) = save2(i)/a(i,1)
          end if

        end do

        do i = nde+1,n
          a(i,1) = 0.0E+00
        end do

      end if

    end if

    do i = 1,nde
      save1(i) = save2(i)/ywt(i)
    end do

    sum = snrm2(nde, save1, 1)
    sum0 = 1.0E+00 / max(1.0, abs(t))
    smax = max(sum0, sum)
    smin = min(sum0, sum)
    sum = smax*sqrt(1.0E+00 + (smin/smax)**2)/sqrt(real(nde))
    h = sign(min(2.0*eps/sum, abs(h)), h)
    yh(1:n,2) = h * save2(1:n)

    if (miter == 2 .or. miter == 5 .or. iswflg == 3) then
      do i = 1,n
        fac(i) = sqrt(uround)
      end do
    end if

  else

    if (miter /= mtrold) then
      mtrold = miter
      rc = 0.0E+00
      convrg = .false.
    end if

    if (mint /= mntold) then
      mntold = mint
      oldl0 = el(1,nq)
      call sdcst (maxord, mint, iswflg,  el, tq)
      rc = rc*el(1,nq)/oldl0
      nwait = nq + 2
    end if

    if (h /= hold) then
      nwait = nq + 2
      rh = h/hold
      call sdscl (hmax, n, nq, rmax,  hold, rc, rh, yh)
    end if

  end if

  return
end
subroutine sdntp(h,k,n,nq,t,tout,yh,y)
!
!*******************************************************************************
!
!! SDNTP interpolates the k-th derivative of y at tout,
!   using the data in the yh array.  if k has a value greater than nq,
!   the nq-th derivative is calculated.
!
  integer n
!
  real factor
  real h
  integer i
  integer j
  integer jj
  integer k
  integer kk
  integer kused
  integer nq
  real r
  real t
  real tout
  real y(*)
  real yh(n,*)
!
  if (k == 0) then

    y(1:n) = yh(1:n,nq+1)

    r = ( tout - t ) / h

    do jj = 1, nq
      j = nq + 1 - jj
      y(1:n) = yh(1:n,j) + r*y(1:n)
    end do

  else

    kused = min(k, nq)
    factor = 1.0E+00
    do kk = 1,kused
      factor = factor*real(nq+1-kk)
    end do

    y(1:n) = factor * yh(1:n,nq+1)

    do jj = kused+1,nq

      j = k + 1 + nq - jj
      factor = 1.0E+00

      do kk = 1,kused
        factor = factor*real(j-kk)
      end do

      y(1:n) = factor * yh(1:n,j) + r * y(1:n)

    end do

    y(1:n) = y(1:n) * h**(-kused)

  end if
end
function sdot ( n, x, incx, y, incy )
!
!*******************************************************************************
!
!! SDOT forms the dot product of two vectors.
!
!
!  Modified:
!
!    02 June 2000
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, real X(*), one of the vectors to be multiplied.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input, real Y(*), one of the vectors to be multiplied.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
!    Output, real SDOT, the dot product of X and Y.
!
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real sdot
  real stemp
  real x(*)
  real y(*)
!
  if ( n <= 0 ) then

    sdot = 0.0E+00

  else if ( incx == 1 .and. incy == 1 ) then

    sdot = dot_product ( x(1:n), y(1:n) )

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    stemp = 0.0E+00
    do i = 1, n
      stemp = stemp + x(ix) * y(iy)
      ix = ix + incx
      iy = iy + incy
    end do

    sdot = stemp

  end if

  return
end
subroutine sdpsc(ksgn,n,nq,yh)
!
!*******************************************************************************
!
!! SDPSC computes the predicted yh values by effectively
!     multiplying the yh array by the pascal triangle matrix when ksgn
!     is +1, and performs the inverse function when ksgn is -1.
!
  integer n
!
  integer i
  integer j
  integer j1
  integer j2
  integer ksgn
  integer nq
  real yh(n,*)
!
  if ( ksgn > 0 ) then

    do j1 = 1, nq
      do j2 = j1, nq
        j = nq - j2 + j1
        do i = 1, n
          yh(i,j) = yh(i,j) + yh(i,j+1)
        end do
      end do
    end do

  else

    do j1 = 1, nq
      do j2 = j1, nq
        j = nq - j2 + j1
        do i = 1, n
          yh(i,j) = yh(i,j) - yh(i,j+1)
        end do
      end do
    end do
  end if

  return
end
subroutine sdpst (el,f,fa,h,impl,jacobn,matdim,miter,ml,mu,n,nde, &
  nq,save2,t,users,y,yh,ywt,uround,nfe,nje,a,dfdy,fac,ier,ipvt, &
  save1,iswflg,bnd,jstate)
!
!*******************************************************************************
!
!! SDPST is called to reevaluate the partials.
!
!
!  if miter is 1, 2, 4, or 5, the matrix
!  p = i - l(0)*h*jacobian is stored in dfdy and subjected to lu
!  decomposition, with the results also stored in dfdy.
!
  integer matdim
  integer n
!
  real a(matdim,*)
  real bl
  real bnd
  real bp
  real br
  real bu
  real dfdy(matdim,*)
  real dfdymx
  real diff
  real dy
  real el(13,12)
  real fac(*)
  real facmin
  real, parameter :: facmax = 0.5E+00
  real factor
  real h
  integer i
  integer i1
  integer i2
  integer i3
  logical ier
  integer iflag
  integer imax
  integer impl
  integer info
  integer ipvt(*)
  integer iswflg
  integer j
  integer j2
  integer jstate
  integer k
  integer miter
  integer ml
  integer mu
  integer mw
  integer nde
  integer nfe
  integer nje
  integer nq
  real save1(*)
  real save2(*)
  real scale
  real snrm2
  real t
  real uround
  real y(*)
  real yh(n,*)
  real yj
  real ys
  real ywt(*)
!
  external f
  external fa
  external jacobn
  external users
!
  nje = nje + 1
  ier = .false.
  if (miter == 1 .or. miter == 2) then
    if (miter == 1) then
      call jacobn (n, t, y, dfdy, matdim, ml, mu)
      if (n == 0) then
        jstate = 8
        return
      end if
      if (iswflg == 3) bnd = snrm2(n*n, dfdy, 1)
      factor = -el(1,nq)*h

      do j = 1,n
        do i = 1,n
          dfdy(i,j) = factor*dfdy(i,j)
        end do
      end do

    else if (miter == 2) then
      br = uround**(.875)
      bl = uround**(.75)
      bu = uround**(.25)
      bp = uround**(-.15)
      facmin = uround**(.78)

      do j = 1,n

        ys = max(abs(ywt(j)), abs(y(j)))
 120        dy = fac(j)*ys

        if (dy == 0.0) then
          if (fac(j) < facmax) then
            fac(j) = min(100.0*fac(j), facmax)
            go to 120
          else
            dy = ys
          end if
        end if

        if (nq == 1) then
          dy = sign(dy, save2(j))
        else
          dy = sign(dy, yh(j,3))
        end if

        dy = (y(j) + dy) - y(j)
        yj = y(j)
        y(j) = y(j) + dy
        call f (n, t, y, save1)

        if (n == 0) then
          jstate = 6
          return
        end if

        y(j) = yj
        factor = -el(1,nq)*h/dy
        dfdy(1:n,j) = ( save1(1:n) - save2(1:n) ) * factor
        diff = abs(save2(1) - save1(1))
        imax = 1

        do i = 2,n
          if (abs(save2(i) - save1(i)) > diff) then
            imax = i
            diff = abs(save2(i) - save1(i))
          end if
        end do
!                                                                 step 2
        if (min(abs(save2(imax)), abs(save1(imax))) > 0.0) then
          scale = max(abs(save2(imax)), abs(save1(imax)))
!                                                                 step 3
          if (diff > bu*scale) then
            fac(j) = max(facmin, fac(j)*0.1)
          else if (br*scale <= diff .and. diff <= bl*scale) then
            fac(j) = min(fac(j)*10.0, facmax)
!                                                                 step 4
          else if (diff < br*scale) then
            fac(j) = min(bp*fac(j), facmax)
          end if
        end if

      end do

      if (iswflg == 3) bnd = snrm2(n*n, dfdy, 1)/(-el(1,nq)*h)
      nfe = nfe + n
    end if

    if (impl == 0) then

      do i = 1,n
        dfdy(i,i) = dfdy(i,i) + 1.0E+00
      end do

    else if (impl == 1) then

      call fa (n, t, y, a, matdim, ml, mu, nde)

      if (n == 0) then
        jstate = 9
        return
      end if

      do j = 1,n
        do i = 1,n
          dfdy(i,j) = dfdy(i,j) + a(i,j)
        end do
      end do

    else if (impl == 2) then

      call fa (n, t, y, a, matdim, ml, mu, nde)

      if (n == 0) then
        jstate = 9
        return
      end if

      do i = 1,nde
        dfdy(i,i) = dfdy(i,i) + a(i,1)
      end do
    end if

    call sgefa (dfdy, matdim, n, ipvt, info)
    if (info /= 0) ier = .true.

  else if (miter == 4 .or. miter == 5) then

    if (miter == 4) then

      call jacobn (n, t, y, dfdy(ml+1,1), matdim, ml, mu)

      if (n == 0) then
        jstate = 8
        return
      end if

      factor = -el(1,nq)*h
      mw = ml + mu + 1

      do j = 1,n
        i1 = max(ml+1, mw+1-j)
        i2 = min(mw+n-j, mw+ml)
        do i = i1,i2
          dfdy(i,j) = factor*dfdy(i,j)
        end do
      end do

    else if (miter == 5) then

      br = uround**(.875)
      bl = uround**(.75)
      bu = uround**(.25)
      bp = uround**(-.15)
      facmin = uround**(.78)
      mw = ml + mu + 1
      j2 = min(mw, n)

      do j = 1,j2

        do k = j,n,mw

          ys = max(abs(ywt(k)), abs(y(k)))

 280      continue

          dy = fac(k)*ys

          if (dy == 0.0) then
            if (fac(k) < facmax) then
              fac(k) = min(100.0*fac(k), facmax)
              go to 280
            else
              dy = ys
            end if
          end if

          if (nq == 1) then
            dy = sign(dy, save2(k))
          else
            dy = sign(dy, yh(k,3))
          end if
          dy = (y(k) + dy) - y(k)
          dfdy(mw,k) = y(k)
          y(k) = y(k) + dy

        end do

        call f (n, t, y, save1)

        if (n == 0) then
          jstate = 6
          return
        end if

        do k = j,n,mw

          y(k) = dfdy(mw,k)
          ys = max(abs(ywt(k)), abs(y(k)))
          dy = fac(k)*ys
          if (dy == 0.0) dy = ys
          if (nq == 1) then
            dy = sign(dy, save2(k))
          else
            dy = sign(dy, yh(k,3))
          end if
          dy = (y(k) + dy) - y(k)
          factor = -el(1,nq)*h/dy
          i1 = max(ml+1, mw+1-k)
          i2 = min(mw+n-k, mw+ml)
          do i = i1,i2
            i3 = k + i - mw
            dfdy(i,k) = factor*(save1(i3) - save2(i3))
          end do

          imax = max(1, k - mu)
          diff = abs(save2(imax) - save1(imax))
          i1 = imax
          i2 = min(k + ml, n)

          do i = i1+1,i2
            if (abs(save2(i) - save1(i)) > diff) then
              imax = i
              diff = abs(save2(i) - save1(i))
            end if
          end do

          if (min(abs(save2(imax)), abs(save1(imax))) >0.0) then

            scale = max(abs(save2(imax)), abs(save1(imax)))

            if (diff > bu*scale) then
              fac(k) = max(facmin, fac(k)*.1)
            else if (br*scale <=diff .and. diff <=bl*scale) then
              fac(k) = min(fac(k)*10.0, facmax)
            else if (diff < br*scale) then
              fac(k) = min(bp*fac(k), facmax)
            end if

          end if

        end do

      end do

      nfe = nfe + j2

    end if

    if (iswflg == 3) then

      dfdymx = 0.0E+00

      do j = 1,n
        i1 = max(ml+1, mw+1-j)
        i2 = min(mw+n-j, mw+ml)
        do i = i1,i2
          dfdymx = max(dfdymx, abs(dfdy(i,j)))
        end do
      end do

      bnd = 0.0E+00
      if (dfdymx /= 0.0) then
        do j = 1,n
          i1 = max(ml+1, mw+1-j)
          i2 = min(mw+n-j, mw+ml)
          do i = i1,i2
                bnd = bnd + (dfdy(i,j)/dfdymx)**2
          end do
        end do
        bnd = dfdymx*sqrt(bnd)/(-el(1,nq)*h)
      end if

    end if

    if (impl == 0) then

      dfdy(mw,1:n) = dfdy(mw,1:n) + 1.0E+00

    else if (impl == 1) then

      call fa (n, t, y, a(ml+1,1), matdim, ml, mu, nde)

      if (n == 0) then
        jstate = 9
        return
      end if

      do j = 1,n
        i1 = max(ml+1, mw+1-j)
        i2 = min(mw+n-j, mw+ml)
        do i = i1,i2
          dfdy(i,j) = dfdy(i,j) + a(i,j)
        end do
      end do

    else if (impl == 2) then

      call fa (n, t, y, a, matdim, ml, mu, nde)

      if (n == 0) then
        jstate = 9
        return
      end if

      do j = 1,nde
        dfdy(mw,j) =  dfdy(mw,j) + a(j,1)
      end do

    end if

    call sgbfa (dfdy, matdim, n, ml, mu, ipvt, info)
    if (info /= 0) ier = .true.

  else if (miter == 3) then

    iflag = 1
    call users ( y, yh(1,2), ywt, save1, save2, t, h, el(1,nq), impl, n, &
      nde, iflag)

    if (n == 0) then
      jstate = 10
      return
    end if

  end if

  return
end
subroutine sdriv1 (n,t,y,tout,mstate,eps,work,lenw)
!
!*******************************************************************************
!
!! SDRIV1 solves ordinary differential equations of the form
!            dy(i)/dt = f(y(i),t), given the initial conditions
!            y(i) = yi.  sdriv1 uses single precision arithmetic.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  i.  choosing the correct routine
!
!     sdriv
!     ddriv
!     cdriv
!           these are the generic names for three packages for solving
!           initial value problems for ordinary differential equations.
!           sdriv uses single precision arithmetic.  ddriv uses double
!           precision arithmetic.  cdriv allows complex-valued
!           differential equations, integrated with respect to a single,
!           real, independent variable.
!
!    as an aid in selecting the proper program, the following is a
!    discussion of the important options or restrictions associated with
!    each program:
!
!      a. sdriv1 should be tried first for those routine problems with
!         no more than 200 differential equations.  internally this
!         routine has two important technical defaults:
!           1. numerical approximation of the jacobian matrix of the
!              right hand side is used.
!           2. the stiff solver option is used.
!         most users of sdriv1 should not have to concern themselves
!         with these details.
!
!      b. sdriv2 should be considered for those problems for which
!         sdriv1 is inadequate (sdriv2 has no explicit restriction on
!         the number of differential equations.)  for example, sdriv1
!         may have difficulty with problems having zero initial
!         conditions and zero derivatives.  in this case sdriv2, with an
!         appropriate value of the parameter ewt, should perform more
!         efficiently.  sdriv2 provides three important additional
!         options:
!           1. the nonstiff equation solver (as well as the stiff
!              solver) is available.
!           2. the root-finding option is available.
!           3. the program can dynamically select either the non-stiff
!              or the stiff methods.
!         internally this routine also defaults to the numerical
!         approximation of the jacobian matrix of the right hand side.
!
!      c. sdriv3 is the most flexible, and hence the most complex, of
!         the programs.  its important additional features include:
!           1. the ability to exploit band structure in the jacobian
!              matrix.
!           2. the ability to solve some implicit differential
!              equations, i.e., those having the form:
!                   a(y,t)*dy/dt = f(y,t).
!           3. the option of integrating in the one step mode.
!           4. the option of allowing the user to provide a routine
!              which computes the analytic jacobian matrix of the right
!              hand side.
!           5. the option of allowing the user to provide a routine
!              which does all the matrix algebra associated with
!              corrections to the solution components.
!
!  ii.  abstract
!
!    the function of sdriv1 is to solve n (200 or fewer) ordinary
!    differential equations of the form dy(i)/dt = f(y(i),t), given the
!    initial conditions y(i) = yi.  sdriv1 is to be called once for each
!    output point.
!
!  iii.  parameters
!
!    the user should use parameter names in the call sequence of sdriv1
!    for those quantities whose value may be altered by sdriv1.  the
!    parameters in the call sequence are:
!
!    n      = (input) the number of differential equations, n <= 200
!
!    t      = the independent variable.  on input for the first call, t
!             is the initial point.  on output, t is the point at which
!             the solution is given.
!
!    y      = the vector of dependent variables.  y is used as input on
!             the first call, to set the initial values.  on output, y
!             is the computed solution vector.  this array y is passed
!             in the call sequence of the user-provided routine f.  thus
!             parameters required by f can be stored in this array in
!             components n+1 and above.  (note: changes by the user to
!             the first n components of this array will take effect only
!             after a restart, i.e., after setting mstate to +1(-1).)
!
!    tout   = (input) the point at which the solution is desired.
!
!    mstate = an integer describing the status of integration.  the user
!             must initialize mstate to +1 or -1.  if mstate is
!             positive, the routine will integrate past tout and
!             interpolate the solution.  this is the most efficient
!             mode.  if mstate is negative, the routine will adjust its
!             internal step to reach tout exactly (useful if a
!             singularity exists beyond tout.)  the meaning of the
!             magnitude of mstate:
!               1  (input) means the first call to the routine.  this
!                  value must be set by the user.  on all subsequent
!                  calls the value of mstate should be tested by the
!                  user.  unless sdriv1 is to be reinitialized, only the
!                  sign of mstate may be changed by the user.  (as a
!                  convenience to the user who may wish to put out the
!                  initial conditions, sdriv1 can be called with
!                  mstate=+1(-1), and tout=t.  in this case the program
!                  will return with mstate unchanged, i.e.,
!                  mstate=+1(-1).)
!               2  (output) means a successful integration.  if a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance tout and call
!                  again.  all other parameters are automatically set.
!               3  (output)(unsuccessful) means the integrator has taken
!                  1000 steps without reaching tout.  the user can
!                  continue the integration by simply calling sdriv1
!                  again.
!               4  (output)(unsuccessful) means too much accuracy has
!                  been requested.  eps has been increased to a value
!                  the program estimates is appropriate.  the user can
!                  continue the integration by simply calling sdriv1
!                  again.
!               5  (output)(unsuccessful) n has been set to zero in
!                subroutine f.  see description of f in section iv.
!
!    eps    = on input, the requested relative accuracy in all solution
!             components.  on output, the adjusted relative accuracy if
!             the input value was too small.  the value of eps should be
!             set as large as is reasonable, because the amount of work
!             done by sdriv1 increases as eps decreases.
!
!    work
!    lenw   = (input)
!             work is an array of lenw real words used
!             internally for temporary storage.  the user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       real work(...)
!             the length of work should be at least n*n + 11*n + 225
!             and lenw should be set to the value used.  the contents of
!             work should not be disturbed between calls to sdriv1.
!
!  long description
!
!  iv.  usage
!
!                   program sample
!                   real alfa, eps, t, tout
!c                                          n is the number of equations
!                   parameter(alfa = 1.0, n = 3,
!                  8          lenw = n*n + 11*n + 225)
!                   real work(lenw), y(n+1)
!c                                                         initial point
!                   t = 0.00001
!c                                                set initial conditions
!                   y(1) = 10.0E+00
!                   y(2) = 0.0E+00
!                   y(3) = 10.0E+00
!c                                                        pass parameter
!                   y(4) = alfa
!                   tout = t
!                   mstate = 1
!                   eps = .001
!              10   call sdriv1 (n, t, y, tout, mstate, eps, work, lenw)
!                   if (mstate > 2) stop
!                   write(*, '(4e12.3)') tout, (y(i), i=1,3)
!                   tout = 10.0*tout
!                   if (tout < 50.0) go to 10
!                 end
!
!    the user must write a subroutine called f to evaluate the right
!    hand side of the differential equations.  it is of the form:
!
!                 subroutine f (n, t, y, ydot)
!                   real alfa, t, y(*), ydot(*)
!                   alfa = y(n+1)
!                   ydot(1) = 1.0E+00 + alfa*(y(2) - y(1)) - y(1)*y(3)
!                   ydot(2) = alfa*(y(1) - y(2)) - y(2)*y(3)
!                   ydot(3) = 1.0E+00 - y(3)*(y(1) + y(2))
!                 end
!
!    this computes ydot = f(y,t), the right hand side of the
!    differential equations.  here y is a vector of length at least n.
!    the actual length of y is determined by the user's declaration in
!    the program which calls sdriv1.  thus the dimensioning of y in f,
!    while required by fortran convention, does not actually allocate
!    any storage.  when this subroutine is called, the first n
!    components of y are intermediate approximations to the solution
!    components.  the user should not alter these values.  here ydot is
!    a vector of length n.  the user should only compute ydot(i) for i
!    from 1 to n.  normally a return from f passes control back to
!    sdriv1.  however, if the user would like to abort the calculation,
!    i.e., return control to the program which calls sdriv1, he should
!    set n to zero.  sdriv1 will signal this by returning a value of
!    mstate equal to +5(-5).  altering the value of n in f has no effect
!    on the value of n in the call sequence of sdriv1.
!
!  v.  other communication to the user
!
!    a. the solver communicates to the user through the parameters
!       above.  in addition it writes diagnostic messages through the
!       standard error handling program xerror.  that program will
!       terminate the user's run if it detects a probable problem setup
!       error, e.g., insufficient storage allocated by the user for the
!       work array.  for further information see section iii-a of the
!       writeup for sdriv3.
!
!    b. the number of evaluations of the right hand side can be found
!       in the work array in the location determined by:
!            lenw - (n + 21) + 4
!
!  references  gear, c. w., "numerical initial value problems in
!                 ordinary differential equations", prentice-hall, 1971.
!
  integer, parameter :: idliw = 21
  integer, parameter :: mxn = 200
  integer n
!
  real eps
  real ewt(1)
  real hmax
  integer i
  integer, parameter :: ierror = 2
  integer ii
  integer, parameter :: impl = 0
  integer iwork(idliw+mxn)
  integer leniw
  integer lenw
  integer lenwcm
  integer lnwchk
  integer, parameter :: mint = 2
  integer, parameter :: miter = 2
  integer ml
  character msg*103
  integer mstate
  real mu
  integer, parameter :: mxord = 5
  integer, parameter :: mxstep = 1000
  integer nde
  integer, parameter :: nroot = 0
  integer nstate
  integer ntask
  real t
  real tout
  real work(*)
  real y(n)
!
  external f
!
  ewt(1) = 1.0E+00

  if ( n > mxn ) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV1 - Fatal error!'
    write ( *, * ) '  The number of equations is too large.'
    write ( *, * ) '  The number of equations N = ', n
    write ( *, * ) '  The maximum is MXN = ', mxn
    stop
  end if

  if (mstate > 0) then
    nstate = mstate
    ntask = 1
  else
    nstate = - mstate
    ntask = 3
  end if

  hmax = 2.0E+00 * abs ( tout - t )
  leniw = n + idliw
  lenwcm = lenw - leniw

  if (lenwcm < (n*n + 10*n + 204)) then
    lnwchk = n*n + 10*n + 204 + leniw
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV1 - Fatal error!'
    write ( *, * ) '  Insufficient work storage.'
    write ( *, * ) '  The given work storage is = ', lenwcm
    write ( *, * ) '  The required work storage is = ', lnwchk
    stop
  end if

  if ( nstate /= 1 ) then
    do i = 1,leniw
      ii = i + lenwcm
      iwork(i) = int(work(ii))
    end do
  end if

  call sdriv3 (n, t, y, f, nstate, tout, ntask, nroot, eps, ewt, &
    ierror, mint, miter, impl, ml, mu, mxord, hmax, work, &
    lenwcm, iwork, leniw, f, f, nde, mxstep, f, f)

  do i = 1, leniw
    ii = lenwcm + i
    work(ii) = real(iwork(i))
  end do

  if ( nstate <= 4 ) then
    mstate = sign(nstate, mstate)
  else if (nstate == 6) then
    mstate = sign(5, mstate)
  end if

  return
end
subroutine sdriv2 (n,t,y,f,tout,mstate,nroot,eps,ewt,mint,work,lenw,iwork, &
  leniw,g)
!
!*******************************************************************************
!
!! SDRIV2 solves n ordinary differential
!            equations of the form dy(i)/dt = f(y(i),t), given the
!            initial conditions y(i) = yi.  the program has options to
!            allow the solution of both stiff and non-stiff differential
!            equations.  sdriv2 uses single precision arithmetic.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!  i.  abstract
!
!    the function of sdriv2 is to solve n ordinary differential
!    equations of the form dy(i)/dt = f(y(i),t), given the initial
!    conditions y(i) = yi.  the program has options to allow the
!    solution of both stiff and non-stiff differential equations.
!    sdriv2 is to be called once for each output point of t.
!
!  ii.  parameters
!
!    the user should use parameter names in the call sequence of sdriv2
!    for those quantities whose value may be altered by sdriv2.  the
!    parameters in the call sequence are:
!
!    n      = (input) the number of differential equations.
!
!    t      = the independent variable.  on input for the first call, t
!             is the initial point.  on output, t is the point at which
!             the solution is given.
!
!    y      = the vector of dependent variables.  y is used as input on
!             the first call, to set the initial values.  on output, y
!             is the computed solution vector.  this array y is passed
!             in the call sequence of the user-provided routines f and
!             g.  thus parameters required by f and g can be stored in
!             this array in components n+1 and above.  (note: changes
!             by the user to the first n components of this array will
!             take effect only after a restart, i.e., after setting
!             mstate to +1(-1).)
!
!    f      = a subroutine supplied by the user.  the name must be
!             declared external in the user's calling program.  this
!           subroutine is of the form:
!                 subroutine f (n, t, y, ydot)
!                   real y(*), ydot(*)
!                     .
!                     .
!                   ydot(1) = ...
!                     .
!                     .
!                   ydot(n) = ...
!                 end (sample)
!             this computes ydot = f(y,t), the right hand side of the
!             differential equations.  here y is a vector of length at
!             least n.  the actual length of y is determined by the
!             user's declaration in the program which calls sdriv2.
!             thus the dimensioning of y in f, while required by fortran
!             convention, does not actually allocate any storage.  when
!             this subroutine is called, the first n components of y are
!             intermediate approximations to the solution components.
!             the user should not alter these values.  here ydot is a
!             vector of length n.  the user should only compute ydot(i)
!             for i from 1 to n.  normally a return from f passes
!             control back to  sdriv2.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls sdriv2, he should set n to zero.
!             sdriv2 will signal this by returning a value of mstate
!             equal to +6(-6).  altering the value of n in f has no
!             effect on the value of n in the call sequence of sdriv2.
!
!    tout   = (input) the point at which the solution is desired.
!
!    mstate = an integer describing the status of integration.  the user
!             must initialize mstate to +1 or -1.  if mstate is
!             positive, the routine will integrate past tout and
!             interpolate the solution.  this is the most efficient
!             mode.  if mstate is negative, the routine will adjust its
!             internal step to reach tout exactly (useful if a
!             singularity exists beyond tout.)  the meaning of the
!             magnitude of mstate:
!               1  (input) means the first call to the routine.  this
!                  value must be set by the user.  on all subsequent
!                  calls the value of mstate should be tested by the
!                  user.  unless sdriv2 is to be reinitialized, only the
!                  sign of mstate may be changed by the user.  (as a
!                  convenience to the user who may wish to put out the
!                  initial conditions, sdriv2 can be called with
!                  mstate=+1(-1), and tout=t.  in this case the program
!                  will return with mstate unchanged, i.e.,
!                  mstate=+1(-1).)
!               2  (output) means a successful integration.  if a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance tout and call
!                  again.  all other parameters are automatically set.
!               3  (output)(unsuccessful) means the integrator has taken
!                  1000 steps without reaching tout.  the user can
!                  continue the integration by simply calling sdriv2
!                  again.  other than an error in problem setup, the
!                  most likely cause for this condition is trying to
!                  integrate a stiff set of equations with the non-stiff
!                  integrator option. (see description of mint below.)
!               4  (output)(unsuccessful) means too much accuracy has
!                  been requested.  eps has been increased to a value
!                  the program estimates is appropriate.  the user can
!                  continue the integration by simply calling sdriv2
!                  again.
!               5  (output) a root was found at a point less than tout.
!                  the user can continue the integration toward tout by
!                  simply calling sdriv2 again.
!               6  (output)(unsuccessful) n has been set to zero in
!                subroutine f.
!               7  (output)(unsuccessful) n has been set to zero in
!                function g.  see description of g below.
!
!    nroot  = (input) the number of equations whose roots are desired.
!             if nroot is zero, the root search is not active.  this
!             option is useful for obtaining output at points which are
!             not known in advance, but depend upon the solution, e.g.,
!             when some solution component takes on a specified value.
!             the root search is carried out using the user-written
!           function g (see description of g below.)  sdriv2 attempts
!             to find the value of t at which one of the equations
!             changes sign.  sdriv2 can find at most one root per
!             equation per internal integration step, and will then
!             return the solution either at tout or at a root, whichever
!             occurs first in the direction of integration.  the index
!             of the equation whose root is being reported is stored in
!             the sixth element of iwork.
!             note: nroot is never altered by this program.
!
!    eps    = on input, the requested relative accuracy in all solution
!             components.  eps = 0 is allowed.  on output, the adjusted
!             relative accuracy if the input value was too small.  the
!             value of eps should be set as large as is reasonable,
!             because the amount of work done by sdriv2 increases as
!             eps decreases.
!
!    ewt    = (input) problem zero, i.e., the smallest physically
!             meaningful value for the solution.  this is used inter-
!             nally to compute an array ywt(i) = max(abs(y(i)), ewt).
!             one step error estimates divided by ywt(i) are kept less
!             than eps.  setting ewt to zero provides pure relative
!             error control.  however, setting ewt smaller than
!             necessary can adversely affect the running time.
!
!    mint   = (input) the integration method flag.
!               mint = 1  means the adams methods, and is used for
!                         non-stiff problems.
!               mint = 2  means the stiff methods of gear (i.e., the
!                         backward differentiation formulas), and is
!                         used for stiff problems.
!               mint = 3  means the program dynamically selects the
!                         adams methods when the problem is non-stiff
!                         and the gear methods when the problem is
!                         stiff.
!             mint may not be changed without restarting, i.e., setting
!             the magnitude of mstate to 1.
!
!    work
!    lenw   = (input)
!             work is an array of lenw real words used
!             internally for temporary storage.  the user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       real work(...)
!             the length of work should be at least
!               16*n + 2*nroot + 204         if mint is 1, or
!               n*n + 10*n + 2*nroot + 204   if mint is 2, or
!               n*n + 17*n + 2*nroot + 204   if mint is 3,
!             and lenw should be set to the value used.  the contents of
!             work should not be disturbed between calls to sdriv2.
!
!    iwork
!    leniw  = (input)
!             iwork is an integer array of length leniw used internally
!             for temporary storage.  the user must allocate space for
!             this array in the calling program by a statement such as
!                       integer iwork(...)
!             the length of iwork should be at least
!               21      if mint is 1, or
!               n+21    if mint is 2 or 3,
!             and leniw should be set to the value used.  the contents
!             of iwork should not be disturbed between calls to sdriv2.
!
!    g      = a real fortran function supplied by the user
!             if nroot is not 0.  in this case, the name must be
!             declared external in the user's calling program.  g is
!             repeatedly called with different values of iroot to
!             obtain the value of each of the nroot equations for which
!             a root is desired.  g is of the form:
!                   real function g (n, t, y, iroot)
!                   real y(*)
!                   go to (10, ...), iroot
!              10   g = ...
!                     .
!                     .
!                 end (sample)
!             here, y is a vector of length at least n, whose first n
!             components are the solution components at the point t.
!             the user should not alter these values.  the actual length
!             of y is determined by the user's declaration in the
!             program which calls sdriv2.  thus the dimensioning of y in
!             g, while required by fortran convention, does not actually
!             allocate any storage.  normally a return from g passes
!             control back to  sdriv2.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls sdriv2, he should set n to zero.
!             sdriv2 will signal this by returning a value of mstate
!             equal to +7(-7).  in this case, the index of the equation
!             being evaluated is stored in the sixth element of iwork.
!             altering the value of n in g has no effect on the value of
!             n in the call sequence of sdriv2.
!
!  long description
!
!  iii.  other communication to the user
!
!    a. the solver communicates to the user through the parameters
!       above.  in addition it writes diagnostic messages through the
!       standard error handling program xerror.  that program will
!       terminate the user's run if it detects a probable problem setup
!       error, e.g., insufficient storage allocated by the user for the
!       work array.  messages are written on the standard error message
!       file.  at installations which have this error handling package
!       the user should determine the standard error handling file from
!       the local documentation.  otherwise the short but serviceable
!       routine, xerror, available with this package, can be used.  that
!       program writes on logical unit 6 to transmit messages.  a
!       complete description of xerror is given in the sandia
!       laboratories report sand78-1189 by r. e. jones.
!
!    b. the first three elements of work and the first five elements of
!       iwork will contain the following statistical data:
!         avgh     the average step size used.
!         hused    the step size last used (successfully).
!         avgord   the average order used.
!         imxerr   the index of the element of the solution vector that
!                  contributed most to the last error test.
!         nqused   the order last used (successfully).
!         nstep    the number of steps taken since last initialization.
!         nfe      the number of evaluations of the right hand side.
!         nje      the number of evaluations of the jacobian matrix.
!
!  iv.  remarks
!
!    a. on any return from sdriv2 all information necessary to continue
!       the calculation is contained in the call sequence parameters,
!       including the work arrays.  thus it is possible to suspend one
!       problem, integrate another, and then return to the first.
!
!    b. if this package is to be used in an overlay situation, the user
!       must declare in the primary overlay the variables in the call
!       sequence to sdriv2.
!
!    c. when the routine g is not required, difficulties associated with
!       an unsatisfied external can be avoided by using the name of the
!       routine which calculates the right hand side of the differential
!       equations in place of g in the call sequence of sdriv2.
!
!  v.  usage
!
!               program sample
!               external f
!               parameter(mint = 1, nroot = 0, n = ...,
!              8          lenw = 16*n + 2*nroot + 204, leniw = 21)
!                                           n is the number of equations
!               real eps, ewt, t, tout, work(lenw), y(n)
!               integer iwork(leniw)
!               open(file='tape6', unit=6, status='new')
!               t = 0.                           initial point
!               y(1:n) = ...                     set initial conditions
!               tout = t
!               ewt = ...
!               mstate = 1
!               eps = ...
!          20   call sdriv2 (n, t, y, f, tout, mstate, nroot, eps, ewt,
!              8             mint, work, lenw, iwork, leniw, f)
!                                          last argument is not the same
!                                          as f if rootfinding is used.
!               if (mstate > 2) stop
!               write(6, 100) tout, (y(i), i=1,n)
!               tout = tout + 1.
!               if (tout <= 10.) go to 20
!          100  format(...)
!             end (sample)
!
!  references  gear, c. w., "numerical initial value problems in
!                 ordinary differential equations", prentice-hall, 1971.
!
  integer n
!
  real eps
  real ewt
  real ewtcom(1)
  real g
  real hmax
  integer ierror
  integer, parameter :: impl = 0
  integer iwork(*)
  integer leniw
  integer lenw
  integer mint
  integer miter
  integer ml
  character msg*81
  integer mstate
  integer mu
  integer mxord
  integer, parameter :: mxstep = 1000
  integer nde
  integer nroot
  integer nstate
  integer ntask
  real t
  real tout
  real work(*)
  real y(n)
!
  external f
  external g
!
  if (mint < 1 .or. mint > 3) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV2 - Fatal error!'
    write ( *, * ) '  Improper value for the integration method.'
    write ( *, * ) '  MINT = ', mint
    write ( *, * ) '  MINT should be between 1 and 3.'
    stop
  end if

  if (mstate >= 0) then
    nstate = mstate
    ntask = 1
  else
    nstate = - mstate
    ntask = 3
  end if

  ewtcom(1) = ewt

  if (ewt /= 0.0) then
    ierror = 3
  else
    ierror = 2
  end if

  if (mint == 1) then
    miter = 0
    mxord = 12
  else if (mint == 2) then
    miter = 2
    mxord = 5
  else if (mint == 3) then
    miter = 2
    mxord = 12
  end if

  hmax = 2.0*abs(tout - t)

  call sdriv3 (n, t, y, f, nstate, tout, ntask, nroot, eps, ewtcom, &
    ierror, mint, miter, impl, ml, mu, mxord, hmax, work, &
    lenw, iwork, leniw, f, f, nde, mxstep, g, f)

  if (mstate >= 0) then
    mstate = nstate
  else
    mstate = - nstate
  end if

  return
end
subroutine sdriv3 (n,t,y,f,nstate,tout,ntask,nroot,eps,ewt,ierror, &
  mint,miter,impl,ml,mu,mxord,hmax,work,lenw,iwork,leniw,jacobn, &
  fa,nde,mxstep,g,users)
!
!*******************************************************************************
!
!! SDRIV3 solves n ordinary differential
!            equations of the form dy(i)/dt = f(y(i),t), given the
!            initial conditions y(i) = yi.  the program has options to
!            allow the solution of both stiff and non-stiff differential
!            equations.  other important options are available.  sdriv3
!            uses single precision arithmetic.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!  i.  abstract
!
!    the primary function of sdriv3 is to solve n ordinary differential
!    equations of the form dy(i)/dt = f(y(i),t), given the initial
!    conditions y(i) = yi.  the program has options to allow the
!    solution of both stiff and non-stiff differential equations.  in
!    addition, sdriv3 may be used to solve:
!      1. the initial value problem, a*dy(i)/dt = f(y(i),t), where a is
!         a non-singular matrix depending on y and t.
!      2. the hybrid differential/algebraic initial value problem,
!         a*dy(i)/dt = f(y(i),t), where a is a vector (whose values may
!         depend upon y and t) some of whose components will be zero
!         corresponding to those equations which are algebraic rather
!         than differential.
!    sdriv3 is to be called once for each output point of t.
!
!  ii.  parameters
!
!    the user should use parameter names in the call sequence of sdriv3
!    for those quantities whose value may be altered by sdriv3.  the
!    parameters in the call sequence are:
!
!    n      = (input) the number of dependent functions whose solution
!             is desired.  n must not be altered during a problem.
!
!    t      = the independent variable.  on input for the first call, t
!             is the initial point.  on output, t is the point at which
!             the solution is given.
!
!    y      = the vector of dependent variables.  y is used as input on
!             the first call, to set the initial values.  on output, y
!             is the computed solution vector.  this array y is passed
!             in the call sequence of the user-provided routines f,
!             jacobn, fa, users, and g.  thus parameters required by
!             those routines can be stored in this array in components
!             n+1 and above.  (note: changes by the user to the first
!             n components of this array will take effect only after a
!             restart, i.e., after setting nstate to 1 .)
!
!    f      = a subroutine supplied by the user.  the name must be
!             declared external in the user's calling program.  this
!           subroutine is of the form:
!                 subroutine f (n, t, y, ydot)
!                   real y(*), ydot(*)
!                     .
!                     .
!                   ydot(1) = ...
!                     .
!                     .
!                   ydot(n) = ...
!                 end (sample)
!             this computes ydot = f(y,t), the right hand side of the
!             differential equations.  here y is a vector of length at
!             least n.  the actual length of y is determined by the
!             user's declaration in the program which calls sdriv3.
!             thus the dimensioning of y in f, while required by fortran
!             convention, does not actually allocate any storage.  when
!             this subroutine is called, the first n components of y are
!             intermediate approximations to the solution components.
!             the user should not alter these values.  here ydot is a
!             vector of length n.  the user should only compute ydot(i)
!             for i from 1 to n.  normally a return from f passes
!             control back to  sdriv3.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls sdriv3, he should set n to zero.
!             sdriv3 will signal this by returning a value of nstate
!             equal to 6 .  altering the value of n in f has no effect
!             on the value of n in the call sequence of sdriv3.
!
!    nstate = an integer describing the status of integration.  the
!             meaning of nstate is as follows:
!               1  (input) means the first call to the routine.  this
!                  value must be set by the user.  on all subsequent
!                  calls the value of nstate should be tested by the
!                  user, but must not be altered.  (as a convenience to
!                  the user who may wish to put out the initial
!                  conditions, sdriv3 can be called with nstate=1, and
!                  tout=t.  in this case the program will return with
!                  nstate unchanged, i.e., nstate=1.)
!               2  (output) means a successful integration.  if a normal
!                  continuation is desired (i.e., a further integration
!                  in the same direction), simply advance tout and call
!                  again.  all other parameters are automatically set.
!               3  (output)(unsuccessful) means the integrator has taken
!                  mxstep steps without reaching tout.  the user can
!                  continue the integration by simply calling sdriv3
!                  again.
!               4  (output)(unsuccessful) means too much accuracy has
!                  been requested.  eps has been increased to a value
!                  the program estimates is appropriate.  the user can
!                  continue the integration by simply calling sdriv3
!                  again.
!               5  (output) a root was found at a point less than tout.
!                  the user can continue the integration toward tout by
!                  simply calling sdriv3 again.
!               6  (output)(unsuccessful) n has been set to zero in
!                subroutine f.
!               7  (output)(unsuccessful) n has been set to zero in
!                function g.  see description of g below.
!               8  (output)(unsuccessful) n has been set to zero in
!                subroutine jacobn.  see description of jacobn below.
!               9  (output)(unsuccessful) n has been set to zero in
!                subroutine fa.  see description of fa below.
!              10  (output)(unsuccessful) n has been set to zero in
!                subroutine users.  see description of users below.
!
!    tout   = (input) the point at which the solution is desired.  the
!             position of tout relative to t on the first call
!             determines the direction of integration.
!
!    ntask  = (input) an index specifying the manner of returning the
!             solution, according to the following:
!               ntask = 1  means sdriv3 will integrate past tout and
!                          interpolate the solution.  this is the most
!                          efficient mode.
!               ntask = 2  means sdriv3 will return the solution after
!                          each internal integration step, or at tout,
!                          whichever comes first.  in the latter case,
!                          the program integrates exactly to tout.
!               ntask = 3  means sdriv3 will adjust its internal step to
!                          reach tout exactly (useful if a singularity
!                          exists beyond tout.)
!
!    nroot  = (input) the number of equations whose roots are desired.
!             if nroot is zero, the root search is not active.  this
!             option is useful for obtaining output at points which are
!             not known in advance, but depend upon the solution, e.g.,
!             when some solution component takes on a specified value.
!             the root search is carried out using the user-written
!             function g (see description of g below.)  sdriv3 attempts
!             to find the value of t at which one of the equations
!             changes sign.  sdriv3 can find at most one root per
!             equation per internal integration step, and will then
!             return the solution either at tout or at a root, whichever
!             occurs first in the direction of integration.  the index
!             of the equation whose root is being reported is stored in
!             the sixth element of iwork.
!             note: nroot is never altered by this program.
!
!    eps    = on input, the requested relative accuracy in all solution
!             components.  eps = 0 is allowed.  on output, the adjusted
!             relative accuracy if the input value was too small.  the
!             value of eps should be set as large as is reasonable,
!             because the amount of work done by sdriv3 increases as eps
!             decreases.
!
!    ewt    = (input) problem zero, i.e., the smallest, nonzero,
!             physically meaningful value for the solution.  (array,
!             possibly of length one.  see following description of
!             ierror.)  setting ewt smaller than necessary can adversely
!             affect the running time.
!
!    ierror = (input) error control indicator.  a value of 3 is
!             suggested for most problems.  other choices and detailed
!             explanations of ewt and ierror are given below for those
!             who may need extra flexibility.
!
!             these last three input quantities eps, ewt and ierror
!             control the accuracy of the computed solution.  ewt and
!             ierror are used internally to compute an array ywt.  one
!             step error estimates divided by ywt(i) are kept less than
!             eps in root mean square norm.
!                 ierror (set by the user) =
!                 1  means ywt(i) = 1. (absolute error control)
!                                   ewt is ignored.
!                 2  means ywt(i) = abs(y(i)),  (relative error control)
!                                   ewt is ignored.
!                 3  means ywt(i) = max(abs(y(i)), ewt(1)).
!                 4  means ywt(i) = max(abs(y(i)), ewt(i)).
!                    this choice is useful when the solution components
!                    have differing scales.
!                 5  means ywt(i) = ewt(i).
!             if ierror is 3, ewt need only be dimensioned one.
!             if ierror is 4 or 5, the user must dimension ewt at least
!             n, and set its values.
!
!    mint   = (input) the integration method indicator.
!               mint = 1  means the adams methods, and is used for
!                         non-stiff problems.
!               mint = 2  means the stiff methods of gear (i.e., the
!                         backward differentiation formulas), and is
!                         used for stiff problems.
!               mint = 3  means the program dynamically selects the
!                         adams methods when the problem is non-stiff
!                         and the gear methods when the problem is
!                         stiff.  when using the adams methods, the
!                         program uses a value of miter=0; when using
!                         the gear methods, the program uses the value
!                         of miter provided by the user.  only a value
!                         of impl = 0 and a value of miter = 1, 2, 4, or
!                         5 is allowed for this option.  the user may
!                         not alter the value of mint or miter without
!                         restarting, i.e., setting nstate to 1.
!
!    miter  = (input) the iteration method indicator.
!               miter = 0  means functional iteration.  this value is
!                          suggested for non-stiff problems.
!               miter = 1  means chord method with analytic jacobian.
!                          in this case, the user supplies subroutine
!                          jacobn (see description below).
!               miter = 2  means chord method with jacobian calculated
!                          internally by finite differences.
!               miter = 3  means chord method with corrections computed
!                          by the user-written routine users (see
!                          description of users below.)  this option
!                          allows all matrix algebra and storage
!                          decisions to be made by the user.  when using
!                          a value of miter = 3, the subroutine fa is
!                          not required, even if impl is not 0.  for
!                          further information on using this option, see
!                          section iv-e below.
!               miter = 4  means the same as miter = 1 but the a and
!                          jacobian matrices are assumed to be banded.
!               miter = 5  means the same as miter = 2 but the a and
!                          jacobian matrices are assumed to be banded.
!
!    impl   = (input) the implicit method indicator.
!               impl = 0 means solving dy(i)/dt = f(y(i),t).
!               impl = 1 means solving a*dy(i)/dt = f(y(i),t),
!                        non-singular a (see description of fa below.)
!                        only mint = 1 or 2, and miter = 1, 2, 3, 4, or
!                        5 are allowed for this option.
!               impl = 2 means solving certain systems of hybrid
!                        differential/algebraic equations (see
!                        description of fa below.)  only mint = 2 and
!                        miter = 1, 2, 3, 4, or 5, are allowed for this
!                        option.
!               the value of impl must not be changed during a problem.
!
!    ml     = (input) the lower half-bandwidth in the case of a banded
!             a or jacobian matrix.  (i.e., maximum(r-c) for nonzero
!             a(r,c).)
!
!    mu     = (input) the upper half-bandwidth in the case of a banded
!             a or jacobian matrix.  (i.e., maximum(c-r).)
!
!    mxord  = (input) the maximum order desired. this is <= 12 for
!             the adams methods and <= 5 for the gear methods.  normal
!             value is 12 and 5, respectively.  if mint is 3, the
!             maximum order used will be min(mxord, 12) when using the
!             adams methods, and min(mxord, 5) when using the gear
!             methods.  mxord must not be altered during a problem.
!
!    hmax   = (input) the maximum magnitude of the step size that will
!             be used for the problem.  this is useful for ensuring that
!             important details are not missed.  if this is not the
!             case, a large value, such as the interval length, is
!             suggested.
!
!    work
!    lenw   = (input)
!             work is an array of lenw real words used
!             internally for temporary storage.  the user must allocate
!             space for this array in the calling program by a statement
!             such as
!                       real work(...)
!             the following table gives the required minimum value for
!             the length of work, depending on the value of impl and
!             miter.  lenw should be set to the value used.  the
!             contents of work should not be disturbed between calls to
!             sdriv3.
!
!      impl =   0                   1                   2
!              ---------------------------------------------------------
! miter =  0   (mxord+4)*n +       not allowed         not allowed
!              2*nroot + 204
!
!         1,2  n*n+(mxord+5)*n     2*n*n+(mxord+5)*n   n*n+(mxord+6)*n
!              + 2*nroot + 204     + 2*nroot + 204     + 2*nroot + 204
!
!          3   (mxord+4)*n +       (mxord+4)*n +       (mxord+4)*n +
!              2*nroot + 204       2*nroot + 204       2*nroot + 204
!
!         4,5  (2*ml+mu)*n +       (4*ml+2*mu)*n +     (2*ml+mu)*n +
!              (mxord+6)*n +       (mxord+7)*n +       (mxord+7)*n +
!              2*nroot + 204       2*nroot + 204       2*nroot + 204
!              ---------------------------------------------------------
!
!    iwork
!    leniw  = (input)
!             iwork is an integer array of length leniw used internally
!             for temporary storage.  the user must allocate space for
!             this array in the calling program by a statement such as
!                       integer iwork(...)
!             the length of iwork should be at least
!               21      if miter is 0 or 3, or
!               n+21    if miter is 1, 2, 4, or 5, or mint is 3,
!             and leniw should be set to the value used.  the contents
!             of iwork should not be disturbed between calls to sdriv3.
!
!    jacobn = a subroutine supplied by the user, if miter is 1 or 4.
!             if this is the case, the name must be declared external in
!             the user's calling program.  given a system of n
!             differential equations, it is meaningful to speak about
!             the partial derivative of the i-th right hand side with
!             respect to the j-th dependent variable.  in general there
!             are n*n such quantities.  often however the equations can
!             be ordered so that the i-th differential equation only
!             involves dependent variables with index near i, e.g., i+1,
!             i-2.  such a system is called banded.  if, for all i, the
!             i-th equation depends on at most the variables
!               y(i-ml), y(i-ml+1), ... , y(i), y(i+1), ... , y(i+mu)
!             then we call ml+mu+1 the bandwith of the system.  in a
!             banded system many of the partial derivatives above are
!             automatically zero.  for the cases miter = 1, 2, 4, and 5,
!             some of these partials are needed.  for the cases
!             miter = 2 and 5 the necessary derivatives are
!             approximated numerically by sdriv3, and we only ask the
!             user to tell sdriv3 the value of ml and mu if the system
!             is banded.  for the cases miter = 1 and 4 the user must
!             derive these partials algebraically and encode them in
!           subroutine jacobn.  by computing these derivatives the
!             user can often save 20-30 per cent of the computing time.
!             usually, however, the accuracy is not much affected and
!             most users will probably forego this option.  the optional
!             user-written subroutine jacobn has the form:
!                 subroutine jacobn (n, t, y, dfdy, matdim, ml, mu)
!                   real y(*), dfdy(matdim,*)
!                     .
!                     .
!                     calculate values of dfdy
!                     .
!                     .
!                 end (sample)
!             here y is a vector of length at least n.  the actual
!             length of y is determined by the user's declaration in the
!             program which calls sdriv3.  thus the dimensioning of y in
!             jacobn, while required by fortran convention, does not
!             actually allocate any storage.  when this subroutine is
!             called, the first n components of y are intermediate
!             approximations to the solution components.  the user
!             should not alter these values.  if the system is not
!             banded (miter=1), the partials of the i-th equation with
!             respect to the j-th dependent function are to be stored in
!             dfdy(i,j).  thus partials of the i-th equation are stored
!             in the i-th row of dfdy.  if the system is banded
!             (miter=4), then the partials of the i-th equation with
!             respect to y(j) are to be stored in dfdy(k,j), where
!             k=i-j+mu+1 .  normally a return from jacobn passes control
!             back to sdriv3.  however, if the user would like to abort
!             the calculation, i.e., return control to the program which
!             calls sdriv3, he should set n to zero.  sdriv3 will signal
!             this by returning a value of nstate equal to +8(-8).
!             altering the value of n in jacobn has no effect on the
!             value of n in the call sequence of sdriv3.
!
!    fa     = a subroutine supplied by the user if impl is 1 or 2, and
!             miter is not 3.  if so, the name must be declared external
!             in the user's calling program.  this subroutine computes
!             the array a, where a*dy(i)/dt = f(y(i),t).
!             there are two cases:
!
!               impl=1.
!               subroutine fa is of the form:
!                 subroutine fa (n, t, y, a, matdim, ml, mu, nde)
!                   real y(*), a(matdim,*)
!                     .
!                     .
!                     calculate all values of a
!                     .
!                     .
!                 end (sample)
!               in this case a is assumed to be a nonsingular matrix,
!               with the same structure as dfdy (see jacobn description
!               above).  programming considerations prevent complete
!               generality.  if miter is 1 or 2, a is assumed to be full
!               and the user must compute and store all values of
!               a(i,j), i,j=1, ... ,n.  if miter is 4 or 5, a is assumed
!               to be banded with lower and upper half bandwidth ml and
!               mu.  the left hand side of the i-th equation is a linear
!               combination of dy(i-ml)/dt, dy(i-ml+1)/dt, ... ,
!               dy(i)/dt, ... , dy(i+mu-1)/dt, dy(i+mu)/dt.  thus in the
!               i-th equation, the coefficient of dy(j)/dt is to be
!               stored in a(k,j), where k=i-j+mu+1.
!               note: the array a will be altered between calls to fa.
!
!               impl=2.
!               subroutine fa is of the form:
!                 subroutine fa (n, t, y, a, matdim, ml, mu, nde)
!                   real y(*), a(*)
!                     .
!                     .
!                     calculate non-zero values of a(1),...,a(nde)
!                     .
!                     .
!                 end (sample)
!               in this case it is assumed that the system is ordered by
!               the user so that the differential equations appear
!               first, and the algebraic equations appear last.  the
!               algebraic equations must be written in the form:
!               0 = f(y(i),t).  when using this option it is up to the
!               user to provide initial values for the y(i) that satisfy
!               the algebraic equations as well as possible.  it is
!               further assumed that a is a vector of length nde.  all
!               of the components of a, which may depend on t, y(i),
!               etc., must be set by the user to non-zero values.
!             here y is a vector of length at least n.  the actual
!             length of y is determined by the user's declaration in the
!             program which calls sdriv3.  thus the dimensioning of y in
!             fa, while required by fortran convention, does not
!             actually allocate any storage.  when this subroutine is
!             called, the first n components of y are intermediate
!             approximations to the solution components.  the user
!             should not alter these values.  fa is always called
!             immediately after calling f, with the same values of t
!             and y.  normally a return from fa passes control back to
!             sdriv3.  however, if the user would like to abort the
!             calculation, i.e., return control to the program which
!             calls sdriv3, he should set n to zero.  sdriv3 will signal
!             this by returning a value of nstate equal to +9(-9).
!             altering the value of n in fa has no effect on the value
!             of n in the call sequence of sdriv3.
!
!    nde    = (input) the number of differential equations.  this is
!             required only for impl = 2, with nde < n.
!
!    mxstep = (input) the maximum number of internal steps allowed on
!             one call to sdriv3.
!
!    g      = a real fortran function supplied by the user
!             if nroot is not 0.  in this case, the name must be
!             declared external in the user's calling program.  g is
!             repeatedly called with different values of iroot to obtain
!             the value of each of the nroot equations for which a root
!             is desired.  g is of the form:
!                   real function g (n, t, y, iroot)
!                   real y(*)
!                   go to (10, ...), iroot
!              10   g = ...
!                     .
!                     .
!                 end (sample)
!             here, y is a vector of length at least n, whose first n
!             components are the solution components at the point t.
!             the user should not alter these values.  the actual length
!             of y is determined by the user's declaration in the
!             program which calls sdriv3.  thus the dimensioning of y in
!             g, while required by fortran convention, does not actually
!             allocate any storage.  normally a return from g passes
!             control back to  sdriv3.  however, if the user would like
!             to abort the calculation, i.e., return control to the
!             program which calls sdriv3, he should set n to zero.
!             sdriv3 will signal this by returning a value of nstate
!             equal to +7(-7).  in this case, the index of the equation
!             being evaluated is stored in the sixth element of iwork.
!             altering the value of n in g has no effect on the value of
!             n in the call sequence of sdriv3.
!
!    users  = a subroutine supplied by the user, if miter is 3.
!             if this is the case, the name must be declared external in
!             the user's calling program.  the routine users is called
!             by sdriv3 when certain linear systems must be solved.  the
!             user may choose any method to form, store and solve these
!             systems in order to obtain the solution result that is
!             returned to sdriv3.  in particular, this allows sparse
!             matrix methods to be used.  the call sequence for this
!             routine is:
!
!              subroutine users (y, yh, ywt, save1, save2, t, h, el,
!               8                  impl, n, nde, iflag)
!                real y(*), yh(*), ywt(*), save1(*),
!               8     save2(*), t, h, el
!
!             the input variable iflag indicates what action is to be
!             taken.subroutine users should perform the following
!             operations, depending on the value of iflag and impl.
!
!               iflag = 0
!                 impl = 0.  users is not called.
!                 impl = 1 or 2.  solve the system a*x = save2,
!                   returning the result in save2.  the array save1 can
!                   be used as a work array.
!
!               iflag = 1
!                 impl = 0.  compute, decompose and store the matrix
!                   (i - h*el*j), where i is the identity matrix and j
!                   is the jacobian matrix of the right hand side.  the
!                   array save1 can be used as a work array.
!                 impl = 1 or 2. compute, decompose and store the matrix
!                   (a - h*el*j).  the array save1 can be used as a work
!                   array.
!
!               iflag = 2
!                 impl = 0.   solve the system
!                     (i - h*el*j)*x = h*save2 - yh - save1,
!                   returning the result in save2.
!                 impl = 1 or 2.  solve the system
!                   (a - h*el*j)*x = h*save2 - a*(yh + save1)
!                   returning the result in save2.
!                 the array save1 should not be altered.
!             normally a return from users passes control back to
!             sdriv3.  however, if the user would like to abort the
!             calculation, i.e., return control to the program which
!             calls sdriv3, he should set n to zero.  sdriv3 will signal
!             this by returning a value of nstate equal to +10(-10).
!             altering the value of n in users has no effect on the
!             value of n in the call sequence of sdriv3.
!
!  long description
!
!  iii.  other communication to the user
!
!    a. the solver communicates to the user through the parameters
!       above.  in addition it writes diagnostic messages through the
!       standard error handling program xerror.  that program will
!       terminate the user's run if it detects a probable problem setup
!       error, e.g., insufficient storage allocated by the user for the
!       work array.  messages are written on the standard error message
!       file.  at installations which have this error handling package
!       the user should determine the standard error handling file from
!       the local documentation.  otherwise the short but serviceable
!       routine, xerror, available with this package, can be used.  that
!       program writes on logical unit 6 to transmit messages.  a
!       complete description of xerror is given in the sandia
!       laboratories report sand78-1189 by r. e. jones.  following is a
!       list of possible errors.  unless otherwise noted, all messages
!       come from sdriv3:
!
!        no.  type         message
!        ---  ----         -------
!         1   fatal        from sdriv2: the integration method flag has
!                          an illegal value.
!         2   warning      the output point is inconsistent with the
!                          value of ntask and t.
!         3   warning      number of steps to reach tout exceeds mxstep.
!         4   recoverable  requested accuracy is too stringent.
!         5   warning      step size is below the roundoff level.
!         6   fatal        eps is less than zero.
!         7   fatal        n is not positive.
!         8   fatal        insufficient work space provided.
!         9   fatal        improper value for nstate, mint, miter and/or
!                          impl.
!        10   fatal        the iwork array is too small.
!        11   fatal        the step size has gone to zero.
!        12   fatal        excessive amount of work.
!        13   fatal        for impl=1 or 2, the matrix a is singular.
!        14   fatal        mxord is not positive.
!        15   fatal        from sdriv1: n is greater than 200.
!        16   fatal        from sdriv1: the work array is too small.
!
!    b. the first three elements of work and the first five elements of
!       iwork will contain the following statistical data:
!         avgh     the average step size used.
!         hused    the step size last used (successfully).
!         avgord   the average order used.
!         imxerr   the index of the element of the solution vector that
!                  contributed most to the last error test.
!         nqused   the order last used (successfully).
!         nstep    the number of steps taken since last initialization.
!         nfe      the number of evaluations of the right hand side.
!         nje      the number of evaluations of the jacobian matrix.
!
!  iv.  remarks
!
!    a. other routines used:
!         sdntp, sdzro, sdstp, sdntl, sdpst, sdcor, sdcst,
!         sdpsc, and sdscl;
!         sgefa, sgesl, sgbfa, sgbsl, and snrm2 (from linpack)
!         r1mach (from the bell laboratories machine constants package)
!         xerror (from the slatec common math library)
!       the last seven routines above, not having been written by the
!       present authors, are not explicitly part of this package.
!
!    b. on any return from sdriv3 all information necessary to continue
!       the calculation is contained in the call sequence parameters,
!       including the work arrays.  thus it is possible to suspend one
!       problem, integrate another, and then return to the first.
!
!    c. if this package is to be used in an overlay situation, the user
!       must declare in the primary overlay the variables in the call
!       sequence to sdriv3.
!
!    d. changing parameters during an integration.
!       the value of nroot, eps, ewt, ierror, mint, miter, or hmax may
!       be altered by the user between calls to sdriv3.  for example, if
!       too much accuracy has been requested (the program returns with
!       nstate = 4 and an increased value of eps) the user may wish to
!       increase eps further.  in general, prudence is necessary when
!       making changes in parameters since such changes are not
!       implemented until the next integration step, which is not
!       necessarily the next call to sdriv3.  this can happen if the
!       program has already integrated to a point which is beyond the
!       new point tout.
!
!    e. as the price for complete control of matrix algebra, the sdriv3
!       users option puts all responsibility for jacobian matrix
!       evaluation on the user.  it is often useful to approximate
!       numerically all or part of the jacobian matrix.  however this
!       must be done carefully.  the fortran sequence below illustrates
!       the method we recommend.  it can be inserted directly into
!       subroutine users to approximate jacobian elements in rows i1
!       to i2 and columns j1 to j2.
!              real dfdy(n,n), epsj, h, r, r1mach,
!             8     save1(n), save2(n), t, uround, y(n), yj, ywt(n)
!              uround = r1mach(4)
!              epsj = sqrt(uround)
!              do j = j1,j2
!                r = epsj*max(abs(ywt(j)), abs(y(j)))
!                if (r == 0.0) r = ywt(j)
!                yj = y(j)
!                y(j) = y(j) + r
!                call f (n, t, y, save1)
!                if (n == 0) return
!                y(j) = yj
!                do i = i1,i2
!                  dfdy(i,j) = (save1(i) - save2(i))/r
!                end do
!              end do
!
!       many problems give rise to structured sparse jacobians, e.g.,
!       block banded.  it is possible to approximate them with fewer
!       function evaluations than the above procedure uses; see curtis,
!       powell and reid, j. inst. maths applics, (1974), vol. 13,
!       pp. 117-119.
!
!    f. when any of the routines jacobn, fa, g, or users, is not
!       required, difficulties associated with unsatisfied externals can
!       be avoided by using the name of the routine which calculates the
!       right hand side of the differential equations in place of the
!       corresponding name in the call sequence of sdriv3.
!
!  references  gear, c. w., "numerical initial value problems in
!                 ordinary differential equations", prentice-hall, 1971.
!
  integer n
!
  real ae
  real big
  logical convrg
  real eps
  real ewt(*)
  real g
  real glast
  real h
  real hmax
  real hsign
  integer i
  integer ia
  integer idfdy
  integer ierror
  integer ifac
  integer iflag
  integer ignow
  integer impl
  integer imxerr
  integer info
  integer iroot
  integer isave1
  integer isave2
  integer itroot
  integer iwork(*)
  integer iywt
  integer j
  integer ja
  integer jaml
  integer jerror
  integer jgnow
  integer jhyp
  integer jroot
  integer jsave2
  integer jstate
  integer jtroot
  integer jyh
  integer jywt
  integer lenchk
  integer leniw
  integer lenw
  integer liwchk
  integer matdim
  integer maxord
  integer mint
  integer miter
  integer ml
  character msg*205
  integer mu
  integer mxord
  integer mxstep
  integer nde
  integer ndecom
  integer npar
  integer nroot
  real, parameter :: nround = 20.0E+00
  integer nstate
  integer nstepl
  integer ntask
  real re
  real r1mach
  real size
  real snrm2
  real sum
  real t
  real tlast
  real tout
  real troot
  real uround
  real work(*)
  real y(*)

  integer, parameter :: iavgh = 1
  integer, parameter :: ihused = 2
  integer, parameter :: iavgrd = 3
  integer, parameter :: iel = 4
  integer, parameter :: ih = 160
  integer, parameter :: ihmax = 161
  integer, parameter :: ihold = 162
  integer, parameter :: ihsign = 163
  integer, parameter :: irc = 164
  integer, parameter :: irmax = 165
  integer, parameter :: it = 166
  integer, parameter :: itout = 167
  integer, parameter :: itq = 168
  integer, parameter :: itrend = 204
  integer, parameter :: iyh = 205
  integer, parameter :: indmxr = 1
  integer, parameter :: inqusd = 2
  integer, parameter :: instep = 3
  integer, parameter :: infe = 4
  integer, parameter :: inje = 5
  integer, parameter :: inroot = 6
  integer, parameter :: icnvrg = 7
  integer, parameter :: ijroot = 8
  integer, parameter :: ijtask = 9
  integer, parameter :: imntld = 10
  integer, parameter :: imtrld = 11
  integer, parameter :: inq = 12
  integer, parameter :: inrtld = 13
  integer, parameter :: indtrt = 14
  integer, parameter :: inwait = 15
  integer, parameter :: imnt = 16
  integer, parameter :: imtrsv = 17
  integer, parameter :: imtr = 18
  integer, parameter :: imxrds = 19
  integer, parameter :: imxord = 20
  integer, parameter :: indprt = 21
  integer, parameter :: indpvt = 22
!
  external f, jacobn, fa, g, users
!
  npar = n
  uround = epsilon ( uround )

  if ( nroot /= 0 ) then
    ae = tiny ( ae )
    re = uround
  end if

  if (eps < 0.0) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Fatal error!'
    write ( *, * ) '  Improper value of EPS.'
    write ( *, * ) '  EPS = ', eps
    write ( *, * ) '  EPS should be nonnegative.'
    stop
  end if

  if (n <= 0) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Fatal error!'
    write ( *, * ) '  Improper value for the number of equations.'
    write ( *, * ) '  N = ', n
    write ( *, * ) '  N should be positive.'
    stop
  end if

  if (mxord <= 0) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Fatal error!'
    write ( *, * ) '  Improper value for the maximum order.'
    write ( *, * ) '  MXORD = ', mxord
    write ( *, * ) '  MXORD should be positive.'
    stop
  end if

  if ((mint < 1 .or. mint > 3) .or. (mint == 3 .and. &
       (miter == 0 .or. miter == 3 .or. impl /= 0)) &
       .or. (miter < 0 .or. miter > 5) .or. &
       (impl /= 0 .and. impl /= 1 .and. impl /= 2) .or. &
       ((impl == 1 .or. impl == 2) .and. miter == 0) .or. &
       (impl == 2 .and. mint == 1) .or. &
       (nstate < 1 .or. nstate > 10)) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Fatal error!'
    write ( *, * ) '  Improper value for some input quantity.'
    write ( *, * ) '  NSTATE/MSTATE/MINT/MITER/IMPL.'
    stop
  end if

  if (miter == 0 .or. miter == 3) then
    liwchk = indpvt - 1
  else if (miter == 1 .or. miter == 2 .or. miter == 4 .or. miter == 5) then
    liwchk = indpvt + n - 1
  end if

  if (leniw < liwchk) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Fatal error!'
    write ( *, * ) '  Insufficient integer storage.'
    write ( *, * ) '  LENIW = ', leniw
    write ( *, * ) '  Required = ', liwchk
    stop
  end if
!
!  Allocate the work array
!  iyh is the index of yh in work.
!
  if (mint == 1 .or. mint == 3) then
    maxord = min(mxord, 12)
  else if (mint == 2) then
    maxord = min(mxord, 5)
  end if
  idfdy = iyh + (maxord + 1)*n
!
!  idfdy is the index of dfdy
!
  if (miter == 0 .or. miter == 3)  then
    iywt = idfdy
  else if (miter == 1 .or. miter == 2)  then
    iywt = idfdy + n*n
  else if (miter == 4 .or. miter == 5)  then
    iywt = idfdy + (2*ml + mu + 1)*n
  end if
!
!  iywt is the index of ywt
!
  isave1 = iywt + n
!                                           isave1 is the index of save1
  isave2 = isave1 + n
!                                           isave2 is the index of save2
  ignow = isave2 + n
!                                             ignow is the index of gnow
  itroot = ignow + nroot
!                                           itroot is the index of troot
  ifac = itroot + nroot
!                                               ifac is the index of fac
  if (miter == 2 .or. miter == 5 .or. mint == 3) then
    ia = ifac + n
  else
    ia = ifac
  end if
!
!  ia is the index of a
!
  if (impl == 0 .or. miter == 3) then
    lenchk = ia - 1
  else if (impl == 1 .and. (miter == 1 .or. miter == 2)) then
    lenchk = ia - 1 + n*n
  else if (impl == 1 .and. (miter == 4 .or. miter == 5)) then
    lenchk = ia - 1 + (2*ml + mu + 1)*n
  else if (impl == 2 .and. miter /= 3) then
    lenchk = ia - 1 + n
  end if

  if (lenw < lenchk) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Fatal error!'
    write ( *, * ) '  Insufficient real storage.'
    write ( *, * ) '  LENW = ', lenw
    write ( *, * ) '  Required = ', lenchk
    stop
  end if

  if (miter == 0 .or. miter == 3) then
    matdim = 1
  else if (miter == 1 .or. miter == 2) then
    matdim = n
  else if (miter == 4 .or. miter == 5) then
    matdim = 2*ml + mu + 1
  end if

  if (impl == 0 .or. impl == 1) then
    ndecom = n
  else if (impl == 2) then
    ndecom = nde
  end if

  if (nstate == 1) then
!                                                  initialize parameters
    if (mint == 1 .or. mint == 3) then
      iwork(imxord) = min(mxord, 12)
    else if (mint == 2) then
      iwork(imxord) = min(mxord, 5)
    end if
    iwork(imxrds) = mxord
    if (mint == 1 .or. mint == 2) then
      iwork(imnt) = mint
      iwork(imtr) = miter
      iwork(imntld) = mint
      iwork(imtrld) = miter
    else if (mint == 3) then
      iwork(imnt) = 1
      iwork(imtr) = 0
      iwork(imntld) = iwork(imnt)
      iwork(imtrld) = iwork(imtr)
      iwork(imtrsv) = miter
    end if

    work(ihmax) = hmax
    h = (tout - t)*(1.0E+00 - 4.0*uround)
    h = sign(min(abs(h), hmax), h)
    work(ih) = h
    hsign = sign(1.0, h)
    work(ihsign) = hsign
    iwork(ijtask) = 0
    work(iavgh) = 0.0E+00
    work(ihused) = 0.0E+00
    work(iavgrd) = 0.0E+00
    iwork(indmxr) = 0
    iwork(inqusd) = 0
    iwork(instep) = 0
    iwork(infe) = 0
    iwork(inje) = 0
    iwork(inroot) = 0
    work(it) = t
    iwork(icnvrg) = 0
    iwork(indprt) = 0
!
!  set initial conditions
!
    do i = 1,n
      jyh = i + iyh - 1
      work(jyh) = y(i)
    end do

    if (t == tout) return
    go to 180
  end if
!
!  on a continuation, check that output points have
!  been or will be overtaken.
!
  if (iwork(icnvrg) == 1) then
    convrg = .true.
  else
    convrg = .false.
  end if
  t = work(it)
  h = work(ih)
  hsign = work(ihsign)
  if (iwork(ijtask) == 0) go to 180
!
!                                   iwork(ijroot) flags unreported
!                                   roots, and is set to the value of
!                                   ntask when a root was last selected.
!                                   it is set to zero when all roots
!                                   have been reported.  iwork(inroot)
!                                   contains the index and work(itout)
!                                   contains the value of the root last
!                                   selected to be reported.
!                                   iwork(inrtld) contains the value of
!                                   nroot and iwork(indtrt) contains
!                                   the value of itroot when the array
!                                   of roots was last calculated.
  if (nroot /= 0) then
    jroot = iwork(ijroot)
    if (jroot > 0) then
!
!  tout has just been reported.
!  if troot <= tout, report troot.
!
      if (nstate /= 5) then

       if (tout*hsign >= work(itout)*hsign) then
          troot = work(itout)
          call sdntp(h, 0, n, iwork(inq), t, troot, work(iyh),  y)
          t = troot
          nstate = 5
          go to 580
        end if
!
!  a root has just been reported.
!  select the next root.
!
      else
        troot = t
        iroot = 0

        do i = 1,iwork(inrtld)

          jtroot = iwork(indtrt) + i - 1
          if (work(jtroot)*hsign <= troot*hsign) then
!
!  check for multiple roots.
!
            if (work(jtroot) == work(itout) .and. i > iwork(inroot)) then
              iroot = i
              troot = work(jtroot)
              go to 60
            end if

            if (work(jtroot)*hsign > work(itout)*hsign) then
              iroot = i
              troot = work(jtroot)
            end if

          end if

        end do

 60     continue

        iwork(inroot) = iroot
        work(itout) = troot
        iwork(ijroot) = ntask

        if (ntask == 1) then
          if (iroot == 0) then
            iwork(ijroot) = 0
          else
            if (tout*hsign >= troot*hsign) then
              call sdntp(h, 0, n, iwork(inq), t, troot,work(iyh),y)
              nstate = 5
              t = troot
              go to 580
            end if
          end if
        else if (ntask == 2 .or. ntask == 3) then
!
!  if there are no more roots, or the
!  user has altered tout to be less
!  than a root, set ijroot to zero.
!
          if (iroot == 0 .or. (tout*hsign < troot*hsign)) then
            iwork(ijroot) = 0
          else
            call sdntp(h, 0, n, iwork(inq), t, troot, work(iyh), y)
            nstate = 5
            t = troot
            go to 580
          end if

        end if
      end if
    end if
  end if

  if (ntask == 1) then
    nstate = 2
    if (t*hsign >= tout*hsign) then
      call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
      t = tout
      go to 580
    end if
  else if (ntask == 2) then
!
!  Check if TOUT has been reset.
!
    if (t*hsign > tout*hsign) then
      write ( *, * ) ' '
      write ( *, * ) 'SDRIV3 - Warning!'
      write ( *, * ) '  The input T was beyond TOUT.'
      write ( *, * ) '  The solution was obtained by interpolation.'
      call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
      t = tout
      nstate = 2
      go to 580
    end if
!
!  Determine if tout has been overtaken
!
    if (abs(tout - t)<=nround*uround*max(abs(t), abs(tout))) then
      t = tout
      nstate = 2
      go to 560
    end if
!
!  if there are no more roots to report, report t.
!
    if (nstate == 5) then
      nstate = 2
      go to 560
    end if

    nstate = 2
!                                                       see if tout will
!                                                       be overtaken.
    if ((t + h)*hsign > tout*hsign) then
      h = tout - t
      if ((t + h)*hsign > tout*hsign) h = h*(1.0E+00 - 4.0*uround)
      work(ih) = h
      if (h == 0.0) go to 670
      iwork(ijtask) = -1
    end if

  else if (ntask == 3) then

    nstate = 2

    if (t*hsign > tout*hsign) then
      write ( *, * ) ' '
      write ( *, * ) 'SDRIV3 - Warning!'
      write ( *, * ) '  The input T was beyond TOUT.'
      write ( *, * ) '  The solution was obtained by interpolation.'
      call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
      t = tout
      go to 580
    end if

    if (abs(tout - t)<=nround*uround*max(abs(t), abs(tout))) then
      t = tout
      go to 560
    end if
    if ((t + h)*hsign > tout*hsign) then
      h = tout - t
      if ((t + h)*hsign > tout*hsign) h = h*(1.0E+00 - 4.0*uround)
      work(ih) = h
      if (h == 0.0) go to 670
      iwork(ijtask) = -1
    end if
  end if
!
!  implement changes in mint, miter, and/or hmax.
!
  if ((mint /= iwork(imntld) .or. miter /= iwork(imtrld)) .and. &
    mint /= 3 .and. iwork(imntld) /= 3) iwork(ijtask) = -1
  if (hmax /= work(ihmax)) then
    h = sign(min(abs(h), hmax), h)
    if (h /= work(ih)) then
      iwork(ijtask) = -1
      work(ih) = h
    end if
    work(ihmax) = hmax
  end if
!
 180  nstepl = iwork(instep)

  do i = 1,n
    jyh = iyh + i - 1
    y(i) = work(jyh)
  end do

  if (nroot /= 0) then
    do i = 1,nroot
      jgnow = ignow + i - 1
      work(jgnow) = g (npar, t, y, i)
      if (npar == 0) then
        iwork(inroot) = i
        nstate = 7
        return
      end if
    end do
  end if

  if (ierror == 1) then
    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = 1.0E+00
    end do
    go to 410
  else if (ierror == 5) then
    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = ewt(i)
    end do
    go to 410
  end if
!                                       reset ywt array.  looping point.
 260  if (ierror == 2) then
    do i = 1,n
      if (y(i) == 0.0) go to 290
      jywt = i + iywt - 1
      work(jywt) = abs(y(i))
    end do
    go to 410
 290    if (iwork(ijtask) == 0) then
      call f (npar, t, y, work(isave2))
      if (npar == 0) then
        nstate = 6
        return
      end if
      iwork(infe) = iwork(infe) + 1
      if (miter == 3 .and. impl /= 0) then
        iflag = 0
        call users(y, work(iyh), work(iywt), work(isave1), &
          work(isave2), t, h, work(iel), impl, npar, ndecom, iflag)
        if (npar == 0) then
          nstate = 10
          return
        end if
      else if (impl == 1) then
        if (miter == 1 .or. miter == 2) then
          call fa (npar, t, y, work(ia), matdim, ml, mu, ndecom)
          if (npar == 0) then
            nstate = 9
            return
          end if
          call sgefa (work(ia), matdim, n, iwork(indpvt), info)
          if (info /= 0) go to 690
          call sgesl(work(ia),matdim,n,iwork(indpvt),work(isave2),0)
        else if (miter == 4 .or. miter == 5) then
          jaml = ia + ml
          call fa (npar, t, y, work(jaml), matdim, ml, mu, ndecom)
          if (npar == 0) then
            nstate = 9
            return
          end if
          call sgbfa (work(ia),matdim,n,ml,mu,iwork(indpvt),info)
          if (info /= 0) go to 690
          call sgbsl (work(ia), matdim, n, ml, mu, iwork(indpvt), &
            work(isave2), 0)

        end if

      else if (impl == 2) then

        call fa (npar, t, y, work(ia), matdim, ml, mu, ndecom)
        if (npar == 0) then
          nstate = 9
          return
        end if

        do i = 1,ndecom
          ja = i + ia - 1
          jsave2 = i + isave2 - 1
          if (work(ja) == 0.0) go to 690
          work(jsave2) = work(jsave2)/work(ja)
        end do

      end if

    end if

    do j = i,n

      jywt = j + iywt - 1

      if (y(j) /= 0.0) then
        work(jywt) = abs(y(j))
      else
        if (iwork(ijtask) == 0) then
          jsave2 = j + isave2 - 1
          work(jywt) = abs(h*work(jsave2))
        else
          jhyp = j + iyh + n - 1
          work(jywt) = abs(work(jhyp))
        end if
      end if

      if (work(jywt) == 0.0) work(jywt) = uround

    end do

  else if (ierror == 3) then

    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = max(ewt(1), abs(y(i)))
    end do

  else if (ierror == 4) then

    do i = 1,n
      jywt = i + iywt - 1
      work(jywt) = max(ewt(i), abs(y(i)))
    end do

  end if

 410  continue

  do i = 1,n
    jywt = i + iywt - 1
    jsave2 = i + isave2 - 1
    work(jsave2) = y(i)/work(jywt)
  end do

  sum = snrm2(n, work(isave2), 1)/sqrt(real(n))

  if (eps < sum*uround) then
    eps = sum*uround*(1.0E+00 + 10.0*uround)
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Warning!'
    write ( *, * ) '  The requested accuracy EPS was not obtainable.'
    write ( *, * ) '  EPS has been increased to ', eps
    nstate = 4
    go to 560
  end if

  if (abs(h) >= uround*abs(t)) then
    iwork(indprt) = 0
  else if (iwork(indprt) == 0) then
    write ( *, * ) ' '
    write ( *, * ) 'SDRIV3 - Warning!'
    write ( *, * ) '  The stepsize is smaller than roundoff.'
    write ( *, * ) '  This may occur when there is an abrupt change'
    write ( *, * ) '  in the right hand side.'
    iwork(indprt) = 1
  end if

  if ( ntask /= 2 ) then
    if ((iwork(instep)-nstepl) > mxstep) then
      write ( *, * ) ' '
      write ( *, * ) 'SDRIV3 - Warning!'
      write ( *, * ) '  Number of steps taken = ', mxstep
      write ( *, * ) '  TOUT not reached.'
      nstate = 3
      go to 560
    end if
  end if

  call sdstp (eps, f, fa, work(ihmax), impl, jacobn, matdim, &
    iwork(imxord), iwork(imnt), iwork(imtr), ml, mu, npar, &
    ndecom, work(iywt), uround, users,  work(iavgh), &
    work(iavgrd), work(ih), work(ihused), iwork(ijtask), &
    iwork(imntld), iwork(imtrld), iwork(infe), iwork(inje), &
    iwork(inqusd), iwork(instep), work(it), y, work(iyh), &
    work(ia), convrg, work(idfdy), work(iel), work(ifac), &
    work(ihold), iwork(indpvt), jstate, iwork(inq), &
    iwork(inwait), work(irc), work(irmax), work(isave1), &
    work(isave2), work(itq), work(itrend), mint, &
    iwork(imtrsv), iwork(imxrds))

  t = work(it)
  h = work(ih)
  go to (470, 670, 680, 690, 690, 660, 660, 660, 660, 660), jstate

 470  iwork(ijtask) = 1
!
!  determine if a root has been overtaken
!
  if (nroot /= 0) then
    iroot = 0

    do i = 1,nroot

      jtroot = itroot + i - 1
      jgnow = ignow + i - 1
      glast = work(jgnow)
      work(jgnow) = g (npar, t, y, i)

      if (npar == 0) then
        iwork(inroot) = i
        nstate = 7
        return
      end if

      if (glast*work(jgnow) > 0.0) then
        work(jtroot) = t + h
      else
        if (work(jgnow) == 0.0) then
          work(jtroot) = t
          iroot = i
        else
          if (glast == 0.0) then
            work(jtroot) = t + h
          else
            if (abs(work(ihused)) >= uround*abs(t)) then
              tlast = t - work(ihused)
              iroot = i
              troot = t
              call sdzro (ae, g, h, npar, iwork(inq), iroot, re, t, &
                work(iyh), uround,  troot, tlast, work(jgnow), glast,  y)
              do j = 1,n
                y(j) = work(iyh + j -1)
              end do
              if (npar == 0) then
                iwork(inroot) = i
                nstate = 7
                return
              end if
              work(jtroot) = troot
            else
              work(jtroot) = t
              iroot = i
            end if
          end if
        end if
      end if

    end do

    if (iroot == 0) then
      iwork(ijroot) = 0
!                                                  select the first root
    else

      iwork(ijroot) = ntask
      iwork(inrtld) = nroot
      iwork(indtrt) = itroot
      troot = t + h

      do i = 1,nroot
        jtroot = itroot + i - 1
        if (work(jtroot)*hsign < troot*hsign) then
          troot = work(jtroot)
          iroot = i
        end if
      end do

      iwork(inroot) = iroot
      work(itout) = troot

      if (troot*hsign <= tout*hsign) then
        call sdntp (h, 0, n, iwork(inq), t, troot, work(iyh),  y)
        nstate = 5
        t = troot
        go to 580
      end if

    end if
  end if
!                               test for ntask condition to be satisfied
  nstate = 2
  if (ntask == 1) then
    if (t*hsign < tout*hsign) go to 260
    call sdntp (h, 0, n, iwork(inq), t, tout, work(iyh),  y)
    t = tout
    go to 580
!                               tout is assumed to have been attained
!                               exactly if t is within twenty roundoff
!                               units of tout, relative to max(tout, t).
  else if (ntask == 2) then
    if (abs(tout - t)<=nround*uround*max(abs(t), abs(tout))) then
      t = tout
    else
      if ((t + h)*hsign > tout*hsign) then
        h = tout - t
        if ((t + h)*hsign>tout*hsign) h = h*(1.0E+00 - 4.0*uround)
        work(ih) = h
        if (h == 0.0) go to 670
        iwork(ijtask) = -1
      end if
    end if
  else if (ntask == 3) then
    if (abs(tout - t)<=nround*uround*max(abs(t), abs(tout))) then
      t = tout
    else
      if ((t + h)*hsign > tout*hsign) then
        h = tout - t
        if ((t + h)*hsign>tout*hsign) h = h*(1.0E+00 - 4.0*uround)
        work(ih) = h
        if (h == 0.0) go to 670
        iwork(ijtask) = -1
      end if
      go to 260
    end if
  end if
!  all returns are made through this
!  section.  imxerr is determined.
 560  continue

  do i = 1,n
    jyh = i + iyh - 1
    y(i) = work(jyh)
  end do

 580  if (convrg) then
    iwork(icnvrg) = 1
  else
    iwork(icnvrg) = 0
  end if
  if (iwork(ijtask) == 0) return
  big = 0.0E+00
  imxerr = 1
  iwork(indmxr) = imxerr
  do i = 1,n
!                                            size = abs(error(i)/ywt(i))
    jywt = i + iywt - 1
    jerror = i + isave1 - 1
    size = abs(work(jerror)/work(jywt))
    if (big < size) then
      big = size
      imxerr = i
      iwork(indmxr) = imxerr
    end if
  end do

  return
!
 660  nstate = jstate
  return

670 continue

  write ( *, * ) ' '
  write ( *, * ) 'SDRIV3 - Fatal error!'
  write ( *, * ) '  The attempted stepsize has been reduced to zero.'
  write ( *, * ) '  The problem setup may be incorrect.'
  stop

680 continue

  write ( *, * ) ' '
  write ( *, * ) 'SDRIV3 - Fatal error!'
  write ( *, * ) '  The stepsize has been reduced about 50 times.'
  write ( *, * ) '  The problem setup may be incorrect.'
  stop

690 continue

  write ( *, * ) ' '
  write ( *, * ) 'SDRIV3 - Fatal error!'
  write ( *, * ) '  Matrix A is singular.'
  stop
end
subroutine sdscl (hmax,n,nq,rmax,h,rc,rh,yh)
!
!*******************************************************************************
!
!! SDSCL rescales the YH array whenever the ODE step size is changed.
!
!
!  SDSCL is a utility routine for the SDRIV family of ODE solvers.
!
  integer n
!
  real h
  real hmax
  integer j
  integer nq
  real r1
  real rc
  real rh
  real rmax
  real yh(n,*)
!
  if ( h < 1.0E+00 ) then
    rh = min ( abs(h)*rh, abs(h)*rmax, hmax ) / abs(h)
  else
    rh = min ( rh, rmax, hmax/abs(h) )
  end if

  r1 = 1.0E+00

  do j = 1, nq
    r1 = r1 * rh
    yh(1:n,j+1) = yh(1:n,j+1) * r1
  end do

  h = h*rh
  rc = rc*rh

  return
end
subroutine sdstp (eps,f,fa,hmax,impl,jacobn,matdim,maxord,mint, &
  miter,ml,mu,n,nde,ywt,uround,users,avgh,avgord,h,hused,jtask, &
  mntold,mtrold,nfe,nje,nqused,nstep,t,y,yh,a,convrg,dfdy,el,fac, &
  hold,ipvt,jstate,nq,nwait,rc,rmax,save1,save2,tq,trend,iswflg, &
  mtrsv,mxrdsv)
!
!*******************************************************************************
!
!! SDSTP performs one step of the integration of an initial value
!  problem for a system of ordinary differential equations.
!  communication with sdstp is done with the following variables:
!
!    yh      an n by maxord+1 array containing the dependent variables
!              and their scaled derivatives.  maxord, the maximum order
!              used, is currently 12 for the adams methods and 5 for the
!              gear methods.  yh(i,j+1) contains the j-th derivative of
!              y(i), scaled by h**j/factorial(j).  only y(i),
!              1 <= i <= n, need be set by the calling program on
!              the first entry.  the yh array should not be altered by
!              the calling program.  when referencing yh as a
!              2-dimensional array, use a column length of n, as this is
!              the value used in sdstp.
!    dfdy    a block of locations used for partial derivatives if miter
!              is not 0.  if miter is 1 or 2 its length must be at least
!              n*n.  if miter is 4 or 5 its length must be at least
!              (2*ml+mu+1)*n.
!    ywt     an array of n locations used in convergence and error tests
!    save1
!    save2   arrays of length n used for temporary storage.
!    ipvt    an integer array of length n used by the linear system
!              solvers for the storage of row interchange information.
!    a       a block of locations used to store the matrix a, when using
!              the implicit method.  if impl is 1, a is a matdim by n
!              array.  if miter is 1 or 2 matdim is n, and if miter is 4
!              or 5 matdim is 2*ml+mu+1.  if impl is 2 its length is n.
!    jtask   an integer used on input.
!              it has the following values and meanings:
!                 == 0  perform the first step.  this value enables
!                         the subroutine to initialize itself.
!                > 0  take a new step continuing from the last.
!                         assumes the last step was successful and
!                         user has not changed any parameters.
!                 < 0  take a new step with a new value of h and/or
!                         mint and/or miter.
!    jstate  a completion code with the following meanings:
!                1  the step was successful.
!                2  a solution could not be obtained with h /= 0.
!                3  a solution was not obtained in mxtry attempts.
!                4  for impl /= 0, the matrix a is singular.
!              on a return with jstate > 1, the values of t and
!              the yh array are as of the beginning of the last
!              step, and h is the last step size attempted.
!
  integer matdim
  integer n
!
  real a(matdim,*)
  real avgh
  real avgord
  real bnd
  logical convrg
  real ctest
  real d
  real denom
  real dfdy(matdim,*)
  real d1
  real el(13,12)
  real eps
  real erdn
  real erup
  real etest
  logical evalfa
  logical evaljc
  real fac(*)
  real h
  real hmax
  real hn
  real hold
  real hs
  real hused
  integer i
  logical, save :: ier = .false.
  integer impl
  integer ipvt(*)
  integer iswflg
  integer iter
  integer j
  integer jstate
  integer jtask
  integer maxord
  integer mint
  integer miter
  integer ml
  integer mntold
  integer mtrold
  integer mtrsv
  integer mu
  integer mxrdsv
  integer nde
  integer nfail
  integer nfe
  integer nje
  integer nq
  integer nqused
  integer nstep
  integer nsv
  integer ntry
  real numer
  integer nwait
  real rc
  real rh
  real rh1
  real rh2
  real rh3
  real rmax
  real save1(*)
  real save2(*)
  real snrm2
  logical switch
  real t
  real told
  real tq(3,12)
  real trend
  real uround
  real y(*)
  real yh(n,*)
  real ywt(*)
  real y0nrm

  real, parameter :: bias1 = 1.30
  real, parameter :: bias2 = 1.20
  real, parameter :: bias3 = 1.40
  integer, parameter :: mxfail = 3
  integer, parameter :: mxiter = 3
  integer, parameter :: mxtry = 50
  real, parameter :: rctest = .30
  real, parameter :: rmfail = 2.0E+00
  real, parameter :: rmnorm = 10.0E+00
  real, parameter :: trshld = 1.0E+00
!
  external f, jacobn, fa, users
!
  nsv = n
  bnd = 0.0E+00
  switch = .false.
  ntry = 0
  told = t
  nfail = 0
  if (jtask <= 0) then
    call sdntl (eps, f, fa, hmax, hold, impl, jtask, matdim, &
      maxord, mint, miter, ml, mu, n, nde, save1, t, &
      uround, users, y, ywt,  h, mntold, mtrold, nfe, rc, &
      yh,  a, convrg, el, fac, ier, ipvt, nq, nwait, rh, &
      rmax, save2, tq, trend, iswflg, jstate)
    if (n == 0) go to 440
    if (h == 0.0) go to 400
    if (ier) go to 420
  end if

 100  ntry = ntry + 1
  if (ntry > mxtry) go to 410
  t = t + h
  call sdpsc (1, n, nq,  yh)
  evaljc = ((abs(rc - 1.0) > rctest) .and. (miter /= 0))
  evalfa = .not. evaljc

 110  iter = 0

  y(1:n) = yh(1:n,1)

  call f (n, t, y, save2)

  if (n == 0) then
    jstate = 6
    go to 430
  end if

  nfe = nfe + 1

  if (evaljc .or. ier) then
    call sdpst (el, f, fa, h, impl, jacobn, matdim, miter, ml, &
      mu, n, nde, nq, save2, t, users, y, yh, ywt, uround, &
      nfe, nje,  a, dfdy, fac, ier, ipvt, save1, iswflg, &
      bnd, jstate)

    if (n == 0) go to 430
    if (ier) go to 160
    convrg = .false.
    rc = 1.0E+00
  end if

  save1(1:n) = 0.0E+00
!                      up to mxiter corrector iterations are taken.
!                      convergence is tested by requiring the r.m.s.
!                      norm of changes to be less than eps.  the sum of
!                      the corrections is accumulated in the vector
!                      save1(i).  it is approximately equal to the l-th
!                      derivative of y multiplied by
!                      h**l/(factorial(l-1)*el(l,nq)), and is thus
!                      proportional to the actual errors to the lowest
!                      power of h present (h**l).  the yh array is not
!                      altered in the correction loop.  the norm of the
!                      iterate difference is stored in d.  if
!                      iter > 0, an estimate of the convergence rate
!                      constant is stored in trend, and this is used in
!                      the convergence test.
!
 130  call sdcor (dfdy, el, fa, h, impl, ipvt, matdim, miter, ml, &
        mu, n, nde, nq, t, users, y, yh, ywt,  evalfa, save1, &
        save2,  a, d, jstate)

    if (n == 0) go to 430

  if (iswflg == 3 .and. mint == 1) then

    if (iter == 0) then

      numer = snrm2(n, save1, 1)

      do i = 1,n
        dfdy(1,i) = save1(i)
      end do

      y0nrm = snrm2(n, yh, 1)

    else

      denom = numer

      do i = 1,n
        dfdy(1,i) = save1(i) - dfdy(1,i)
      end do

      numer = snrm2(n, dfdy, matdim)

      if (el(1,nq)*numer <= 100.0*uround*y0nrm) then
        if (rmax == rmfail) then
          switch = .true.
          go to 170
        end if
      end if

      do i = 1,n
        dfdy(1,i) = save1(i)
      end do

      if (denom /= 0.0) then
        bnd = max(bnd, numer/(denom*abs(h)*el(1,nq)))
      end if

    end if
  end if

  if (iter > 0) then
    trend = max(.9*trend, d/d1)
  end if

  d1 = d
  ctest = min(2.0*trend, 1.0)*d
  if (ctest <= eps) go to 170
  iter = iter + 1
  if (iter < mxiter) then
    do i = 1,n
      y(i) = yh(i,1) + el(1,nq)*save1(i)
    end do
    call f (n, t, y, save2)
    if (n == 0) then
      jstate = 6
      go to 430
    end if
    nfe = nfe + 1
    go to 130
  end if
!                     the corrector iteration failed to converge in
!                     mxiter tries.  if partials are involved but are
!                     not up to date, they are reevaluated for the next
!                     try.  otherwise the yh array is retracted to its
!                     values before prediction, and h is reduced, if
!                     possible.  if not, a no-convergence exit is taken.
  if (convrg) then
    evaljc = .true.
    evalfa = .false.
    go to 110
  end if

 160  t = told
  call sdpsc (-1, n, nq,  yh)
  nwait = nq + 2
  if (jtask /= 0 .and. jtask /= 2) rmax = rmfail

  if (iter == 0) then
    rh = .3
  else
    rh = .9*(eps/ctest)**(.2)
  end if

  if (rh*h == 0.0) go to 400
  call sdscl (hmax, n, nq, rmax,  h, rc, rh, yh)
  go to 100
!                          the corrector has converged.  convrg is set
!                          to .true. if partial derivatives were used,
!                          to indicate that they may need updating on
!                          subsequent steps.  the error test is made.
 170  convrg = (miter /= 0)

  do i = 1,nde
    save2(i) = save1(i)/ywt(i)
  end do

  etest = snrm2(nde, save2, 1)/(tq(2,nq)*sqrt(real(nde)))
!
!                           the error test failed.  nfail keeps track of
!                           multiple failures.  restore t and the yh
!                           array to their previous values, and prepare
!                           to try the step again.  compute the optimum
!                           step size for this or one lower order.
  if (etest > eps) then
    t = told
    call sdpsc (-1, n, nq,  yh)
    nfail = nfail + 1
    if (nfail < mxfail) then
      if (jtask /= 0 .and. jtask /= 2) rmax = rmfail
      rh2 = 1.0/(bias2*(etest/eps)**(1.0/real(nq+1)))
      if (nq > 1) then
        do i = 1,nde
          save2(i) = yh(i,nq+1)/ywt(i)
        end do
        erdn = snrm2(nde, save2, 1)/(tq(1,nq)*sqrt(real(nde)))
        rh1 = 1.0/max(1.0, bias1*(erdn/eps)**(1.0/real(nq)))
        if (rh2 < rh1) then
          nq = nq - 1
          rc = rc*el(1,nq)/el(1,nq+1)
          rh = rh1
        else
          rh = rh2
        end if
      else
        rh = rh2
      end if
      nwait = nq + 2
      if (rh*h == 0.0) go to 400
      call sdscl (hmax, n, nq, rmax,  h, rc, rh, yh)
      go to 100
    end if
!                control reaches this section if the error test has
!                failed mxfail or more times.  it is assumed that the
!                derivatives that have accumulated in the yh array have
!                errors of the wrong order.  hence the first derivative
!                is recomputed, the order is set to 1, and the step is
!                retried.
    nfail = 0
    jtask = 2
    do i = 1,n
      y(i) = yh(i,1)
    end do

    call sdntl (eps, f, fa, hmax, hold, impl, jtask, matdim, &
      maxord, mint, miter, ml, mu, n, nde, save1, t, &
      uround, users, y, ywt,  h, mntold, mtrold, nfe, rc, &
      yh,  a, convrg, el, fac, ier, ipvt, nq, nwait, rh, &
      rmax, save2, tq, trend, iswflg, jstate)
    rmax = rmnorm
    if (n == 0) go to 440
    if (h == 0.0) go to 400
    if (ier) go to 420
    go to 100
  end if
!                          after a successful step, update the yh array.
  nstep = nstep + 1
  hused = h
  nqused = nq
  avgh = (real(nstep-1)*avgh + h)/real(nstep)
  avgord = (real(nstep-1)*avgord + real(nq))/real(nstep)

  do j = 1,nq+1
    do i = 1,n
      yh(i,j) = yh(i,j) + el(j,nq)*save1(i)
    end do
  end do

  y(1:n) = yh(1:n,1)
!                                          if iswflg is 3, consider
!                                          changing integration methods.
  if (iswflg == 3) then
    if (bnd /= 0.0) then
      if (mint == 1 .and. nq <= 5) then
        hn = abs(h)/max(uround, (etest/eps)**(1.0/real(nq+1)))
        hn = min(hn, 1.0/(2.0*el(1,nq)*bnd))
        hs = abs(h)/max(uround,(etest/(eps*el(nq+1,1)))**(1.0/real(nq+1)))
        if (hs > 1.2*hn) then
          mint = 2
          mntold = mint
          miter = mtrsv
          mtrold = miter
          maxord = min(mxrdsv, 5)
          rc = 0.0E+00
          rmax = rmnorm
          trend = 1.0E+00
          call sdcst (maxord, mint, iswflg, el, tq)
          nwait = nq + 2
        end if
      else if (mint == 2) then
        hs = abs(h)/max(uround, (etest/eps)**(1.0/real(nq+1)))
        hn = abs(h)/max(uround,(etest*el(nq+1,1)/eps)**(1.0/real(nq+1)))
        hn = min(hn, 1.0/(2.0*el(1,nq)*bnd))
        if (hn >= hs) then
          mint = 1
          mntold = mint
          miter = 0
          mtrold = miter
          maxord = min(mxrdsv, 12)
          rmax = rmnorm
          trend = 1.0E+00
          convrg = .false.
          call sdcst (maxord, mint, iswflg, el, tq)
          nwait = nq + 2
        end if
      end if
    end if
  end if

  if (switch) then
    mint = 2
    mntold = mint
    miter = mtrsv
    mtrold = miter
    maxord = min(mxrdsv, 5)
    nq = min(nq, maxord)
    rc = 0.0E+00
    rmax = rmnorm
    trend = 1.0E+00
    call sdcst (maxord, mint, iswflg, el, tq)
    nwait = nq + 2
  end if
!                           consider changing h if nwait = 1.  otherwise
!                           decrease nwait by 1.  if nwait is then 1 and
!                           nq<maxord, then save1 is saved for use in
!                           a possible order increase on the next step.
!
  if (jtask == 0 .or. jtask == 2) then

    rh = 1.0/max(uround, bias2*(etest/eps)**(1.0/real(nq+1)))
    if (rh>trshld) call sdscl (hmax, n, nq, rmax, h, rc, rh, yh)

  else if (nwait > 1) then

    nwait = nwait - 1

    if (nwait == 1 .and. nq < maxord) then
      do i = 1,nde
        yh(i,maxord+1) = save1(i)
      end do
    end if
!             if a change in h is considered, an increase or decrease in
!             order by one is considered also.  a change in h is made
!             only if it is by a factor of at least trshld.  factors
!             rh1, rh2, and rh3 are computed, by which h could be
!             multiplied at order nq - 1, order nq, or order nq + 1,
!             respectively.  the largest of these is determined and the
!             new order chosen accordingly.  if the order is to be
!             increased, we compute one additional scaled derivative.
!             if there is a change of order, reset nq and the
!             coefficients.  in any case h is reset according to rh and
!             the yh array is rescaled.
  else
    if (nq == 1) then
      rh1 = 0.0E+00
    else
      do i = 1,nde
        save2(i) = yh(i,nq+1)/ywt(i)
      end do
      erdn = snrm2(nde, save2, 1)/(tq(1,nq)*sqrt(real(nde)))
      rh1 = 1.0/max(uround, bias1*(erdn/eps)**(1.0/real(nq)))
    end if

    rh2 = 1.0/max(uround, bias2*(etest/eps)**(1.0/real(nq+1)))

    if (nq == maxord) then
      rh3 = 0.0E+00
    else
      do i = 1,nde
        save2(i) = (save1(i) - yh(i,maxord+1))/ywt(i)
      end do
      erup = snrm2(nde, save2, 1)/(tq(3,nq)*sqrt(real(nde)))
      rh3 = 1.0/max(uround, bias3*(erup/eps)**(1.0/real(nq+2)))
    end if

    if (rh1 > rh2 .and. rh1 >= rh3) then
      rh = rh1
      if (rh <= trshld) go to 380
      nq = nq - 1
      rc = rc*el(1,nq)/el(1,nq+1)
    else if (rh2 >= rh1 .and. rh2 >= rh3) then
      rh = rh2
      if (rh <= trshld) go to 380
    else
      rh = rh3
      if (rh <= trshld) go to 380
      do i = 1,n
        yh(i,nq+2) = save1(i)*el(nq+1,nq)/real(nq+1)
      end do
      nq = nq + 1
      rc = rc*el(1,nq)/el(1,nq-1)
    end if

    if (iswflg == 3 .and. mint == 1) then
      if (bnd/=0.0) rh = min(rh, 1.0/(2.0*el(1,nq)*bnd*abs(h)))
    end if

    call sdscl (hmax, n, nq, rmax,  h, rc, rh, yh)
    rmax = rmnorm
 380    nwait = nq + 2
  end if
!               all returns are made through this section.  h is saved
!               in hold to allow the caller to change h on the next step
  jstate = 1
  hold = h
  return

 400  jstate = 2
  hold = h
  do i = 1,n
    y(i) = yh(i,1)
  end do
  return

 410  jstate = 3
  hold = h
  return

 420  jstate = 4
  hold = h
  return

 430  t = told
  call sdpsc (-1, nsv, nq,  yh)
  do i = 1,nsv
    y(i) = yh(i,1)
  end do
 440  hold = h
  return
end
subroutine sdzro (ae,f,h,n,nq,iroot,re,t,yh,uround,b,c,fb,fc,y)
!
!*******************************************************************************
!
!! SDZRO is a special purpose version of zeroin, modified for use with
!     the sdriv1 package.
!
!     sandia mathematical program library
!     mathematical computing services division 5422
!     sandia laboratories
!     p. o. box 5800
!     albuquerque, new mexico  87115
!     control data 6600 version 4.5, 1 november 1971
!
!  Discussion:
!
!    zeroin searches for a zero of a function f(n, t, y, iroot)
!    between the given values b and c until the width of the
!    interval (b, c) has collapsed to within a tolerance specified
!    by the stopping criterion, abs(b - c) <= 2.*(rw*abs(b) + ae).
!
!  Parameters:
!
!        f     - name of the external function, which returns a
!                real result.  this name must be in an
!                external statement in the calling program.
!        b     - one end of the interval (b, c).  the value returned for
!                b usually is the better approximation to a zero of f.
!        c     - the other end of the interval (b, c).
!        re    - relative error used for rw in the stopping criterion.
!                if the requested re is less than machine precision,
!                then rw is set to approximately machine precision.
!        ae    - absolute error used in the stopping criterion.  if the
!                given interval (b, c) contains the origin, then a
!                nonzero value should be chosen for ae.
!
!     references
!       1.  l f shampine and h a watts, zeroin, a root-solving routine,
!           sc-tm-70-631, sept 1970.
!       2.  t j dekker, finding a zero by means of successive linear
!           interpolation, "constructive aspects of the fundamental
!           theorem of algebra", edited by b dejon and p henrici, 1969.
!
  integer n
!
  real a
  real acbs
  real acmb
  real ae
  real b
  real c
  real cmb
  real er
  real f
  real fa
  real fb
  real fc
  real h
  integer ic
  integer iroot
  integer kount
  integer nq
  real p
  real q
  real re
  real rw
  real t
  real tol
  real uround
  real y(*)
  real yh(n,*)
!
  external f
!
  er = 4.0*uround
  rw = max(re, er)
  ic = 0
  acbs = abs(b - c)
  a = c
  fa = fc
  kount = 0
!                                                    perform interchange
 10 continue

  if (abs(fc) < abs(fb)) then
    a = b
    fa = fb
    b = c
    fb = fc
    c = a
    fc = fa
  end if

  cmb = 0.5*(c - b)
  acmb = abs(cmb)
  tol = rw*abs(b) + ae
!                                                test stopping criterion
  if (acmb <= tol) return
  if (kount > 50) return
!                                    calculate new iterate implicitly as
!                                    b + p/q, where we arrange p >= 0.
!                         the implicit form is used to prevent overflow.
  p = (b - a)*fb
  q = fa - fb
  if (p < 0.0) then
    p = -p
    q = -q
  end if
!                          update a and check for satisfactory reduction
!                          in the size of our bounding interval.
  a = b
  fa = fb
  ic = ic + 1
  if (ic >= 4) then
    if (8.0*acmb >= acbs) then
!                                                                 bisect
      b = 0.5*(c + b)
      go to 20
    end if
    ic = 0
  end if
  acbs = acmb
!                                            test for too small a change
  if (p <= abs(q)*tol) then
!                                                 increment by tolerance
    b = b + sign(tol, cmb)
!                                               root ought to be between
!                                               b and (c + b)/2.
  else if (p < cmb*q) then
!                                                            interpolate
    b = b + p/q
  else
!                                                                 bisect
    b = 0.5*(c + b)
  end if
!
!  have completed computation for new iterate b.
!
 20   continue

  call sdntp (h, 0, n, nq, t, b, yh,  y)
  fb = f(n, b, y, iroot)
  if (n == 0) return
  if (fb == 0.0) return
  kount = kount + 1
!
!  Decide whether next step is interpolation or extrapolation
!
  if (sign(1.0, fb) == sign(1.0, fc)) then
    c = a
    fc = fa
  end if
  go to 10
end
subroutine secfac(nr,n,x,g,a,xpls,gpls,epsm,itncnt,rnf, &
  iagflg,noupdt,s,y,u,w)
!
!*******************************************************************************
!
!! SECFAC updates the hessian by the bfgs factored method.
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> old iterate, x[k-1]
! g(n)         --> gradient or approximate at old iterate
! a(n,n)      <--> on entry: cholesky decomposition of hessian in
!                    lower part and diagonal.
!                  on exit:  updated cholesky decomposition of hessian
!                    in lower triangular part and diagonal
! xpls(n)      --> new iterate, x[k]
! gpls(n)      --> gradient or approximate at new iterate
! epsm         --> machine epsilon
! itncnt       --> iteration count
! rnf          --> relative noise in optimization function fcn
! iagflg       --> =1 if analytic gradient supplied, =0 itherwise
! noupdt      <--> boolean: no update yet
!                  [retain value between successive calls]
! s(n)         --> workspace
! y(n)         --> workspace
! u(n)         --> workspace
! w(n)         --> workspace
!
  integer n
  integer nr
!
  real a(nr,n)
  real alp
  real den1
  real den2
  real epsm
  real g(n)
  real gpls(n)
  integer i
  integer iagflg
  integer itncnt
  integer j
  logical noupdt
  real reltol
  real rnf
  real s(n)
  logical skpupd
  real snorm2
  real snrm2
  real u(n)
  real w(n)
  real x(n)
  real xpls(n)
  real y(n)
  real ynrm2
!
  if ( itncnt==1) noupdt=.true.

  s(1:n) = xpls(1:n) - x(1:n)
  y(1:n) = gpls(1:n) - g(1:n)

  den1 = dot_product ( s, y )

  snorm2 = snrm2(n,s,1)

  ynrm2=snrm2(n,y,1)

  if ( den1<sqrt(epsm)*snorm2*ynrm2) then
    return
  end if

  call mvmltu(nr,n,a,s,u)

  den2 = dot_product ( u, u )
!
!  l <-- sqrt(den1/den2)*l
!
  alp=sqrt(den1/den2)

  if ( noupdt ) then

    do j=1,n
      u(j)=alp*u(j)
      do i=j,n
        a(i,j)=alp*a(i,j)
      end do
    end do

    noupdt=.false.
    den2=den1
    alp=1.0E+00

  end if

  skpupd=.true.
!
!  w = l(l+)s = hs
!
  call mvmltl(nr,n,a,u,w)
  i=1

  if ( iagflg == 0 ) then
    reltol=sqrt(rnf)
  else
    reltol=rnf
  end if

60  continue

  if ( i <= n .and. skpupd ) then

    if ( abs(y(i)-w(i)) < reltol * max(abs(g(i)),abs(gpls(i))) ) then
      i = i + 1
    else
      skpupd=.false.
    end if
    go to 60
  end if

  if ( skpupd ) then
    return
  end if
!
!  w=y-alp*l(l+)s
!
  w(1:n)=y(1:n)-alp*w(1:n)
!
!  alp=1/sqrt(den1*den2)
!
  alp=alp/den1
!
!  u=(l+)/sqrt(den1*den2) = (l+)s/sqrt((y+)s * (s+)l(l+)s)
!
  u(1:n)=alp*u(1:n)
!
!  copy l into upper triangular part.  zero l.
!
  do i=2,n
    do j=1,i-1
      a(j,i)=a(i,j)
      a(i,j)=0.0E+00
    end do
  end do
!
!  find q, (l+) such that  q(l+) = (l+) + u(w+)
!
  call qrupdt(nr,n,a,u,w)
!
!  upper triangular part and diagonal of a now contain updated
!  cholesky decomposition of hessian.  copy back to lower triangular part.
!
  do i=2,n
    do j=1,i-1
      a(i,j)=a(j,i)
    end do
  end do

  return
end
subroutine secunf(nr,n,x,g,a,udiag,xpls,gpls,epsm,itncnt, &
  rnf,iagflg,noupdt,s,y,t)
!
!*******************************************************************************
!
!! SECUNF updates hessian by the bfgs unfactored method
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> old iterate, x[k-1]
! g(n)         --> gradient or approximate at old iterate
! a(n,n)      <--> on entry: approximate hessian at old iterate
!                    in upper triangular part (and udiag)
!                  on exit:  updated approx hessian at new iterate
!                    in lower triangular part and diagonal
!                  [lower triangular part of symmetric matrix]
! udiag        --> on entry: diagonal of hessian
! xpls(n)      --> new iterate, x[k]
! gpls(n)      --> gradient or approximate at new iterate
! epsm         --> machine epsilon
! itncnt       --> iteration count
! rnf          --> relative noise in optimization function fcn
! iagflg       --> =1 if analytic gradient supplied, =0 otherwise
! noupdt      <--> boolean: no update yet
!                  [retain value between successive calls]
! s(n)         --> workspace
! y(n)         --> workspace
! t(n)         --> workspace
!
  integer n
  integer nr
!
  real a(nr,n)
  real den1
  real den2
  real epsm
  real g(n)
  real gam
  real gpls(n)
  integer i
  integer iagflg
  integer itncnt
  integer j
  logical noupdt
  real rnf
  real s(n)
  logical skpupd
  real snorm2
  real snrm2
  real t(n)
  real tol
  real udiag(n)
  real x(n)
  real xpls(n)
  real y(n)
  real ynrm2
!
! copy hessian in upper triangular part and udiag to
! lower triangular part and diagonal
!
  do j = 1, n
    a(j,j) = udiag(j)
    do i = j+1, n
      a(i,j) = a(j,i)
    end do
  end do

  if ( itncnt==1) noupdt=.true.

  s(1:n) = xpls(1:n) - x(1:n)
  y(1:n) = gpls(1:n) - g(1:n)

  den1 = dot_product ( s, y )

  snorm2=snrm2(n,s,1)

  ynrm2=snrm2(n,y,1)

  if ( den1<sqrt(epsm)*snorm2*ynrm2) then
    return
  end if

  call mvmlts(nr,n,a,s,t)

  den2 = dot_product ( s, t )

  if ( noupdt ) then
!
!         h <-- [(s+)y/(s+)hs]h
!
    gam=den1/den2
    den2=gam*den2
    do j=1,n
      t(j)=gam*t(j)
      do i=j,n
        a(i,j)=gam*a(i,j)
      end do
    end do
    noupdt=.false.

  end if

  skpupd=.true.
!
!  check update condition on row i
!
  do i=1,n

    tol=rnf*max(abs(g(i)),abs(gpls(i)))
    if ( iagflg==0) tol=tol/sqrt(rnf)

    if ( abs(y(i)-t(i)) >= tol) then
      skpupd=.false.
      go to 70
    end if

  end do

70 continue

  if ( skpupd) then
    return
  end if
!
!  bfgs update
!
  do j=1,n
    do i=j,n
      a(i,j)=a(i,j)+y(i)*y(j)/den1-t(i)*t(j)/den2
    end do
  end do

  return
end
subroutine sgbfa ( abd, lda, n, ml, mu, ipvt, info )
!
!*******************************************************************************
!
!! SGBFA factors a real band matrix by elimination.
!
!
!  Discussion:
!
!     sgbfa is usually called by sgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!  Parameters:
!
!     on entry
!
!        abd     real(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be >= 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 <= ml < n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 <= mu < n .
!                more efficient if  ml <= mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) == 0.0E+00 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgbsl will divide by zero if
!                     called.  use  rcond  in sgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do j = 1, n
!                      i1 = max ( 1, j-mu )
!                      i2 = min ( n, j+ml )
!                      do i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                      end do
!                   end do
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
  integer lda
  integer n
!
  real abd(lda,n)
  integer i
  integer i0
  integer info
  integer ipvt(n)
  integer isamax
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
  m = ml + mu + 1
  info = 0
!
!  Zero initial fill-in columns.
!
  j0 = mu + 2
  j1 = min ( n, m ) - 1

  do jz = j0, j1
     i0 = m + 1 - jz
     do i = i0, ml
        abd(i,jz) = 0.0E+00
     end do
  end do

  jz = j1
  ju = 0
!
!  Gaussian elimination with partial pivoting.
!
  do k = 1, n-1
!
!  Zero next fill-in column.
!
     jz = jz + 1
     if ( jz <= n ) then
       abd(1:ml,jz) = 0.0E+00
     end if
!
!  Find L = pivot index.
!
     lm = min ( ml, n-k )
     l = isamax ( lm+1, abd(m,k), 1 ) + m - 1
     ipvt(k) = l + k - m
!
!  Zero pivot implies this column already triangularized.
!
     if ( abd(l,k) == 0.0E+00 ) then

       info = k
!
!  Interchange if necessary.
!
     else

        if ( l /= m ) then
           t = abd(l,k)
           abd(l,k) = abd(m,k)
           abd(m,k) = t
        end if
!
!  Compute multipliers.
!
        t = -1.0E+00 / abd(m,k)
        call sscal ( lm, t, abd(m+1,k), 1 )
!
!  Row elimination with column indexing.
!
        ju = min ( max ( ju, mu+ipvt(k) ), n )
        mm = m

        do j = k+1, ju
           l = l - 1
           mm = mm - 1
           t = abd(l,j)
           if ( l /= mm ) then
              abd(l,j) = abd(mm,j)
              abd(mm,j) = t
           end if
           call saxpy ( lm, t, abd(m+1,k), 1, abd(mm+1,j), 1 )
        end do

     end if

  end do

  ipvt(n) = n

  if ( abd(m,n) == 0.0E+00 ) then
    info = n
  end if

  return
end
subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)
!
!*******************************************************************************
!
!! SGBSL solves a banded linear system factored by SGBCO or SGBFA.
!
!
!  Description:
!
!     sgbsl solves the real band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by sbgco or sgbfa.
!
!     on entry
!
!        abd     real(lda, n)
!                the output from sbgco or sgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from sbgco or sgbfa.
!
!        b       real(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if sbgco has set rcond > 0.0E+00
!        or sgbfa has set info == 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sbgco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do j = 1, p
!              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!           end do
!
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!  references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
!                 *linpack users  guide*, siam, 1979.
!
  integer lda
  integer n
!
  real abd(lda,*)
  real b(*)
  integer ipvt(*)
  integer job
  integer k
  integer kb
  integer l
  integer la
  integer lb
  integer lm
  integer m
  integer ml
  integer mu
  real sdot
  real t
!
  m = mu + ml + 1

  if ( job == 0 ) then

    if ( ml /= 0 ) then

      do k = 1, n-1
        lm = min(ml,n-k)
        l = ipvt(k)
        t = b(l)
        if (l /= k) then
          b(l) = b(k)
          b(k) = t
        end if
        call saxpy(lm,t,abd(m+1,k),1,b(k+1),1)
      end do

    end if

    do kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/abd(m,k)
      lm = min(k,m) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      call saxpy(lm,t,abd(la,k),1,b(lb),1)
    end do

  else

     do k = 1, n
       lm = min(k,m) - 1
       la = m - lm
       lb = k - lm
       t = sdot(lm,abd(la,k),1,b(lb),1)
       b(k) = (b(k) - t)/abd(m,k)
     end do

     if ( ml /= 0) then

       do kb = 1, n-1
         k = n - kb
         lm = min(ml,n-k)
         b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
         l = ipvt(k)
         if (l /= k) then
           t = b(l)
           b(l) = b(k)
           b(k) = t
         end if
       end do

     end if

  end if

  return
end
subroutine sgeco(a,lda,n,ipvt,rcond,z)
!
!*******************************************************************************
!
!! SGECO factors a real matrix by gaussian elimination
!     and estimates the condition of the matrix.
!
!     if  rcond  is not needed, sgefa is slightly faster.
!     to solve  a*x = b , follow sgeco by sgesl.
!     to compute  inverse(a)*c , follow sgeco by sgesl.
!     to compute  determinant(a) , follow sgeco by sgedi.
!     to compute  inverse(a) , follow sgeco by sgedi.
!
!     on entry
!
!        a       real(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u , where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        rcond   real
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0E+00 + rcond == 1.0E+00
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z       real(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     linpack.  this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!  references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
!                 *linpack users  guide*, siam, 1979.
!
  integer lda
  integer n
!
  integer ipvt(n)
  real a(lda,*)
  real anorm
  real ek
  integer info
  integer j
  integer k
  integer kb
  integer l
  real rcond
  real s
  real sasum
  real sdot
  real sm
  real t
  real wk
  real wkm
  real ynorm
  real z(*)
!
!  compute 1-norm of a
!
  anorm = 0.0E+00
  do j = 1, n
     anorm = max (anorm,sasum(n,a(1,j),1))
  end do
!
!  factor
!
  call sgefa(a,lda,n,ipvt,info)
!
!     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
!     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
!     trans(a)  is the transpose of a .  the components of  e  are
!     chosen to cause maximum local growth in the elements of w  where
!     trans(u)*w = e .  the vectors are frequently rescaled to avoid
!     overflow.
!
!  solve trans(u)*w = e
!
  ek = 1.0E+00
  z(1:n) = 0.0E+00

  do k = 1, n

     if (z(k) /= 0.0) ek = sign(ek,-z(k))

     if (abs(ek-z(k)) > abs(a(k,k))) then
        s = abs(a(k,k))/abs(ek-z(k))
        call sscal(n,s,z,1)
        ek = s*ek
     end if

     wk = ek - z(k)
     wkm = -ek - z(k)
     s = abs(wk)
     sm = abs(wkm)

     if (a(k,k) /= 0.0) then
       wk = wk/a(k,k)
       wkm = wkm/a(k,k)
     else
       wk = 1.0E+00
       wkm = 1.0E+00
     end if

     if ( k+1 <= n) then

       do j = k+1, n
         sm = sm + abs(z(j)+wkm*a(k,j))
         z(j) = z(j) + wk*a(k,j)
         s = s + abs(z(j))
       end do

       if ( s < sm ) then
           t = wkm - wk
           wk = wkm
           do j = k+1, n
              z(j) = z(j) + t*a(k,j)
           end do
       end if

     end if

     z(k) = wk

  end do

  s = 1.0/sasum(n,z,1)
  call sscal(n,s,z,1)
!
!  solve trans(l)*y = w
!
  do kb = 1, n

     k = n + 1 - kb
     if (k < n) z(k) = z(k) + sdot(n-k,a(k+1,k),1,z(k+1),1)
     if ( abs(z(k)) > 1.0) then
        s = 1.0/abs(z(k))
        call sscal(n,s,z,1)
     end if
     l = ipvt(k)
     t = z(l)
     z(l) = z(k)
     z(k) = t

  end do

  s = 1.0/sasum(n,z,1)
  call sscal(n,s,z,1)

  ynorm = 1.0E+00
!
!  solve l*v = y
!
  do k = 1, n
     l = ipvt(k)
     t = z(l)
     z(l) = z(k)
     z(k) = t
     if (k < n) call saxpy(n-k,t,a(k+1,k),1,z(k+1),1)
     if (abs(z(k)) > 1.0) then
        s = 1.0/abs(z(k))
        call sscal(n,s,z,1)
        ynorm = s*ynorm
     end if
  end do

  s = 1.0/sasum(n,z,1)
  call sscal(n,s,z,1)
  ynorm = s*ynorm
!
!  solve  u*z = v
!
  do kb = 1, n

     k = n + 1 - kb

     if (abs(z(k)) > abs(a(k,k))) then
        s = abs(a(k,k))/abs(z(k))
        call sscal(n,s,z,1)
        ynorm = s*ynorm
     end if

     if (a(k,k) /= 0.0) z(k) = z(k)/a(k,k)
     if (a(k,k) == 0.0) z(k) = 1.0E+00
     t = -z(k)
     call saxpy(k-1,t,a(1,k),1,z(1),1)

  end do
!
!  make znorm = 1.0E+00
!
  s = 1.0/sasum(n,z,1)
  call sscal(n,s,z,1)
  ynorm = s*ynorm

  if (anorm /= 0.0) rcond = ynorm/anorm
  if (anorm == 0.0) rcond = 0.0E+00
  return
end
subroutine sgefa(a,lda,n,ipvt,info)
!
!*******************************************************************************
!
!! SGEFA factors a real matrix by gaussian elimination.
!
!
!     sgefa is usually called by sgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
!
!     on entry
!
!        a       real(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u , where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) == 0.0E+00 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgesl or sgedi will divide by zero
!                     if called.  use  rcond  in sgeco for a reliable
!                     indication of singularity.
!
!  references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
!                 *linpack users  guide*, siam, 1979.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer info
  integer ipvt(*)
  integer isamax
  integer j
  integer k
  integer l
  real t
!
  info = 0

  do k = 1, n-1
!
!  find l = pivot index
!
     l = isamax(n-k+1,a(k,k),1) + k - 1
     ipvt(k) = l
!
!  zero pivot implies this column already triangularized
!
     if (a(l,k) == 0.0E+00 ) then

       info = k
!
!  interchange if necessary
!
    else

        if ( l /= k) then
           t = a(l,k)
           a(l,k) = a(k,k)
           a(k,k) = t
        end if
!
!  compute multipliers
!
        t = -1.0/a(k,k)
        call sscal(n-k,t,a(k+1,k),1)
!
!  row elimination with column indexing
!
        do j = k+1, n
           t = a(l,j)
           if (l /= k) then
             a(l,j) = a(k,j)
             a(k,j) = t
           end if
           call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
        end do

     end if

  end do

  ipvt(n) = n
  if (a(n,n) == 0.0) info = n
  return
end
subroutine sgefs(a,lda,n,v,itask,ind,work,iwork,rcond)
!
!*******************************************************************************
!
!! SGEFS solves a general nxn system of single
!    precision linear equations using linpack subroutines sgeco
!    and sgesl.  that is, if a is an nxn real matrix and if x
!    and b are real n-vectors, then sgefs solves the equation
!
!                          a*x=b.
!
!    the matrix a is first factored into upper and lower tri-
!    angular matrices u and l using partial pivoting.  these
!    factors and the pivoting information are used to find the
!    solution vector x.  an approximate condition number is
!    calculated to provide a rough estimate of the number of
!    digits of accuracy in the computed solution.
!
!    if the equation a*x=b is to be solved for more than one vector
!    b, the factoring of a does not need to be performed again and
!    the option to only solve (itask == 2) will be faster for
!    the succeeding solutions.  in this case, the contents of a,
!    lda, n and iwork must not have been altered by the user follow-
!    ing factorization (itask=1).  ind will not be changed by sgefs
!    in this case.  other settings of itask are used to solve linear
!    systems involving the transpose of a.
!
!  Parameters:
!
!    a      real(lda,n)
!             on entry, the doubly subscripted array with dimension
!               (lda,n) which contains the coefficient matrix.
!             on return, an upper triangular matrix u and the
!               multipliers necessary to construct a matrix l
!               so that a=l*u.
!    lda    integer
!             the leading dimension of the array a.  lda must be great-
!             er than or equal to n.  (terminal error message ind=-1)
!    n      integer
!             the order of the matrix a.  the first n elements of
!             the array a are the elements of the first column of
!             the  matrix a.  n must be greater than or equal to 1.
!             (terminal error message ind=-2)
!    v      real(n)
!             on entry, the singly subscripted array(vector) of di-
!               mension n which contains the right hand side b of a
!               system of simultaneous linear equations a*x=b.
!             on return, v contains the solution vector, x .
!    itask  integer
!             if itask=1, the matrix a is factored and then the
!               linear equation is solved.
!             if itask=2, the equation is solved using the existing
!               factored matrix a and iwork.
!             if itask=3, the matrix is factored and a'x=b is solved
!             if itask=4, the transposed equation is solved using the
!               existing factored matrix a and iwork.
!             if itask < 1 or itask > 4, then the terminal error
!               message ind=-3 is printed.
!    ind    integer
!             gt. 0  ind is a rough estimate of the number of digits
!                     of accuracy in the solution, x.
!             lt. 0  see error message corresponding to ind below.
!    work   real(n)
!             a singly subscripted array of dimension at least n.
!    iwork  integer(n)
!             a singly subscripted array of dimension at least n.
!    rcond  real
!             estimate of 1.0/cond(a)
!
!  error messages printed
!
!    ind=-1  fatal   n is greater than lda.
!    ind=-2  fatal   n is less than 1.
!    ind=-3  fatal   itask is less than 1 or greater than 4.
!    ind=-4  fatal   the matrix a is computationally singular.
!                         a solution has not been computed.
!    ind=-10 warning    the solution has no apparent significance.
!                         the solution may be inaccurate or the matrix
!                         a may be poorly scaled.
!
!  references  subroutine sgefs was developed by group c-3, los alamos
!                 scientific laboratory, los alamos, nm 87545.
!                 the linpack subroutines used by sgefs are described in
!                 detail in the *linpack users guide* published by
!                 the society for industrial and applied mathematics
!                 (siam) dated 1979.
!
  integer lda
  integer n
!
  real a(lda,n)
  integer ind
  integer itask
  integer iwork(*)
  integer job
  real rcond
  real v(*)
  real work(*)
!
  if ( lda < n ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGEFS - Fatal error!'
    write ( *, * ) '  LDA < N.'
    return
  end if

  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGEFS - Fatal error!'
    write ( *, * ) '  N <= 0.'
    return
  end if

  if ( itask < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGEFS - Fatal error!'
    write ( *, * ) '  ITASK < 1.'
    return
  end if

  if ( itask > 4 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SGEFS - Fatal error!'
    write ( *, * ) '  ITASK > 4.'
    return
  end if
!
!  Factor the matrix.
!
  if ( itask == 1 .or. itask == 3 ) then

    call sgeco ( a, lda, n, iwork, rcond, work )
!
!  Check for computational singularity.
!
    if ( rcond == 0.0E+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SGEFS - Error!'
      write ( *, * ) '  The matrix is compuationally singular.'
      return
    end if
!
!  Estimate the number of significant digits.
!
    ind = - int ( log10 ( epsilon ( rcond ) / rcond ) )

    if (ind <= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SGEFS - Warning!'
      write ( *, * ) '  Solution may have no significant digits.'
    end if

  end if

  if ( itask <= 2 ) then
    job = 0
  else
    job = 1
  end if

  call sgesl ( a, lda, n, iwork, v, job )

  return
end
subroutine sgesl(a,lda,n,ipvt,b,job)
!
!*******************************************************************************
!
!! SGESL solves a linear system factored by SGEFA or SGECO.
!
!
!     on entry
!
!        a       real(lda, n)
!                the output from sgeco or sgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from sgeco or sgefa.
!
!        b       real(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically, this indicates singularity,
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if sgeco has set rcond > 0.0E+00
!        or sgefa has set info == 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do j = 1, p
!              call sgesl(a,lda,n,ipvt,c(1,j),0)
!           end do
!
  integer lda
  integer n
!
  real a(lda,n)
  real b(n)
  integer ipvt(*)
  integer job
  integer k
  integer kb
  integer l
  real sdot
  real t
!
  if (job == 0) then
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
     do k = 1, n-1
        l = ipvt(k)
        t = b(l)
        if (l /= k) then
           b(l) = b(k)
           b(k) = t
        end if
        call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
      end do
!
!  solve  u*x = y
!
     do k = n, 1, -1
        b(k) = b(k)/a(k,k)
        t = -b(k)
        call saxpy(k-1,t,a(1,k),1,b(1),1)
     end do

  else
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
     do k = 1, n
       t = sdot(k-1,a(1,k),1,b(1),1)
       b(k) = (b(k) - t)/a(k,k)
     end do
!
!  Solve trans(l)*x = y
!
     do kb = 1, n-1
        k = n - kb
        b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
        l = ipvt(k)
        if (l /= k) then
           t = b(l)
           b(l) = b(k)
           b(k) = t
        end if
     end do

  end if

  return
end
subroutine sinqb ( n, x, wsave )
!
!*******************************************************************************
!
!! SINQB computes the fast sine transform of quarter wave data.
!
!
!  Discussion:
!
!    SINQB computes a sequence from its representation in terms of a sine
!    series with odd wave numbers.
!
!    SINQF is the unnormalized inverse of SINQB since a call of SINQB
!    followed by a call of SINQF will multiply the input sequence X by 4*N.
!
!    The array WSAVE must be initialized by calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!
!        4 * X_in(K) * sin ( ( 2 * K - 1 ) * I * PI / ( 2 * N ) )
!
!  Modified:
!
!    12 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling SINQI.  A different WSAVE array must be used
!    for each different value of N.
!
  integer n
!
  integer k
  real wsave(3*n+15)
  real x(n)
!
  if ( n < 1 ) then
    return
  end if

  if ( n == 1 ) then
    x(1) = 4.0E+00 * x(1)
    return
  end if

  x(2:n:2) = -x(2:n:2)

  call cosqb ( n, x, wsave )
!
!  Reverse the X vector.
!
  call rvec_reverse ( n, x )

  return
end
subroutine sinqf ( n, x, wsave )
!
!*******************************************************************************
!
!! SINQF computes the fast sine transform of quarter wave data.
!
!
!  Discussion:
!
!    SINQF computes the coefficients in a sine series representation with
!    only odd wave numbers.
!
!    SINQB is the unnormalized inverse of SINQF since a call of SINQF
!    followed by a call of SINQB will multiply the input sequence X by 4*N.
!
!    The array WSAVE, which is used by SINQF, must be initialized by
!    calling SINQI.
!
!    The transform is defined by:
!
!      X_out(I) = (-1)**(I-1) * X_in(N) + sum ( 1 <= K <= N-1 )
!        2 * X_in(K) * sin ( ( 2 * I - 1 ) * K * PI / ( 2 * N ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.  The
!    method is more efficient when N is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE(3*N+15), a work array.  The WSAVE array must be
!    initialized by calling SINQI.  A different WSAVE array must be used
!    for each different value of N.
!
  integer n
!
  integer k
  real wsave(3*n+15)
  real x(n)
!
  if ( n <= 1 ) then
    return
  end if
!
!  Reverse the X vector.
!
  call rvec_reverse ( n, x )

  call cosqf ( n, x, wsave )

  x(2:n:2) = -x(2:n:2)

  return
end
subroutine sinqi ( n, wsave )
!
!*******************************************************************************
!
!! SINQI initializes WSAVE, used in SINQF and SINQB.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the array to be transformed.
!
!    Output, real WSAVE(3*N+15), contains data, dependent on the value
!    of N, which is necessary for the SINQF or SINQB routines.
!
  integer n
!
  real wsave(3*n+15)
!
  call cosqi ( n, wsave )

  return
end
subroutine sint ( n, x, wsave )
!
!*******************************************************************************
!
!! SINT computes the discrete Fourier sine transform of an odd sequence.
!
!
!  Discussion:
!
!    SINT is the unnormalized inverse of itself since a call of SINT
!    followed by another call of SINT will multiply the input sequence
!    X by 2*(N+1).
!
!    The array WSAVE must be initialized by calling SINTI.
!
!    The transform is defined by:
!
!      X_out(I) = sum ( 1 <= K <= N )
!        2 * X_in(K) * sin ( K * I * PI / ( N + 1 ) )
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is the product of small primes.
!
!    Input/output, real X(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WSAVE((5*N+30)/2), a work array.  The WSAVE array must be
!    initialized by calling SINTI.  A different WSAVE array must be used
!    for each different value of N.
!
  integer n
!
  integer iw1
  integer iw2
  integer iw3
  real wsave((5*n+30)/2)
  real x(n)
!
  write ( *, * ) 'DEBUG: entered SINT'
  iw1 = n / 2 + 1
  iw2 = iw1 + n + 1
  iw3 = iw2 + n + 1

  call sint1 ( n, x, wsave(1), wsave(iw1), wsave(iw2), wsave(iw3) )
  write ( *, * ) 'DEBUG: leaving SINT'
  return
end
subroutine sint1 ( n, war, was, xh, x, ifac )
!
!*******************************************************************************
!
!! SINT1 is a lower level routine used by SINT.
!
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!
!    Input/output, real WAR(N).
!    On input, the sequence to be transformed.
!    On output, the transformed sequence.
!
!    Input, real WAS(N/2).
!
!    Input, real XH(N).
!
!    Input, real X(N).
!
!    Input, integer IFAC(15).
!    IFAC(1) = N, the number that was factored.
!    IFAC(2) = NF, the number of factors.
!    IFAC(3:2+NF), the factors.
!
  integer n
!
  integer i
  integer ifac(15)
  integer k
  integer kc
  integer ns2
  real, parameter :: sqrt3 = 1.73205080756888E+00
  real t1
  real t2
  real war(n)
  real was(n/2)
  real x(n)
  real xh(n)
  real xhold
!
  write ( *, * ) 'DEBUG: entered SINT1'
  xh(1:n) = war(1:n)
  war(1:n) = x(1:n)

  if ( n <= 1 ) then
    xh(1) = 2.0E+00 * xh(1)
    return
  end if

  if ( n == 2 ) then
    xhold = sqrt3 * ( xh(1) + xh(2) )
    xh(2) = sqrt3 * ( xh(1) - xh(2) )
    xh(1) = xhold
    return
  end if

  ns2 = n / 2
  x(1) = 0.0E+00

  do k = 1, n/2
    t1 = xh(k) - xh(n+1-k)
    t2 = was(k) * ( xh(k) + xh(n+1-k) )
    x(k+1) = t1 + t2
    x(n+2-k) = t2 - t1
  end do

  if ( mod ( n, 2 ) /= 0 ) then
    x(n/2+2) = 4.0E+00 * xh(n/2+1)
  end if
  write ( *, * ) 'DEBUG: calling RFFTF1'
  call rfftf1 ( n+1, x, xh, war, ifac )
  write ( *, * ) 'DEBUG: back from RFFTF1'
  xh(1) = 0.5E+00 * x(1)
  do i = 3, n, 2
    xh(i-1) = -x(i)
    xh(i) = xh(i-2) + x(i-1)
  end do

  if ( mod ( n, 2 ) == 0 ) then
    xh(n) = -x(n+1)
  end if

  x(1:n) = war(1:n)
  war(1:n) = xh(1:n)

  write ( *, * ) 'DEBUG: leaving SINT1'
  return
end
subroutine sinti ( n, wsave )
!
!*******************************************************************************
!
!! SINTI initializes WSAVE, used in SINT.
!
!
!  Discussion:
!
!    The prime factorization of N together with a tabulation of the
!    trigonometric functions are computed and stored in WSAVE.
!
!  Modified:
!
!    09 March 2001
!
!  Author:
!
!    Paul Swarztrauber,
!    National Center for Atmospheric Research
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!    P N Swarztrauber,
!    Vectorizing the FFT's,
!    in Parallel Computations (G. Rodrigue, editor),
!    Academic Press, 1982, pages 51-83.
!
!    B L Buzbee,
!    The SLATEC Common Math Library,
!    in Sources and Development of Mathematical Software (W. Cowell, editor),
!    Prentice Hall, 1984, pages 302-318.
!
!  Parameters:
!
!    Input, integer N, the length of the sequence to be transformed.
!    The method is most efficient when N+1 is a product of small primes.
!
!    Output, real WSAVE((5*N+30)/2), contains data, dependent on the value
!    of N, which is necessary for the SINT routine.
!
  integer n
!
  real dt
  integer k
  real pimach
  real wsave((5*n+30)/2)
!
  if ( n <= 1 ) then
    return
  end if

  dt = pimach() / real ( n + 1 )

  do k = 1, n/2
    wsave(k) = 2.0E+00 * sin ( real ( k ) * dt )
  end do

  call rffti ( n+1, wsave((n/2)+1) )

  return
end
subroutine sndofd(nr,n,xpls,fcn,fpls,a,sx,rnoise,stepsz,anbr)
!
!*******************************************************************************
!
!! SNDOFD approximates a Hessian with a second order finite difference.
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! xpls(n)      --> new iterate:   x[k]
! fcn          --> name of subroutine to evaluate function
! fpls         --> function value at new iterate, f(xpls)
! a(n,n)      <--  finite difference approximation to hessian
!                  only lower triangular matrix and diagonal
!                  are returned
! sx(n)        --> diagonal scaling matrix for x
! rnoise       --> relative noise in fname [f(x)]
! stepsz(n)    --> workspace (stepsize in i-th component direction)
! anbr(n)      --> workspace (neighbor in i-th direction)
!
  integer n
  integer nr
!
  real a(nr,n)
  real anbr(n)
  real fhat
  real fpls
  integer i
  integer j
  real ov3
  real rnoise
  real stepsz(n)
  real sx(n)
  real xpls(n)
  real xtmpi
  real xtmpj
!
  external fcn
!
! find i-th stepsize and evaluate neighbor in direction
! of i-th unit vector.
!
  ov3 = 1.0E+00 / 3.0E+00

  do i = 1, n
    stepsz(i)=rnoise**ov3 * max(abs(xpls(i)),1.0E+00 / sx(i) )
    xtmpi=xpls(i)
    xpls(i)=xtmpi+stepsz(i)
    call fcn(n,xpls,anbr(i))
    xpls(i)=xtmpi
  end do
!
!  calculate column i of a
!
  do i = 1, n

    xtmpi=xpls(i)
    xpls(i)=xtmpi+2.0*stepsz(i)
    call fcn(n,xpls,fhat)
    a(i,i)=((fpls-anbr(i))+(fhat-anbr(i)))/(stepsz(i)*stepsz(i))
!
! calculate sub-diagonal elements of column
!
    if ( i /= n ) then

      xpls(i)=xtmpi+stepsz(i)

      do j=i+1,n
        xtmpj=xpls(j)
        xpls(j)=xtmpj+stepsz(j)
        call fcn(n,xpls,fhat)
        a(j,i)=((fpls-anbr(i))+(fhat-anbr(j)))/(stepsz(i)*stepsz(j))
        xpls(j)=xtmpj
      end do

    end if

    xpls(i) = xtmpi

  end do

  return
end
function snrm2 ( n, x, incx )
!
!*******************************************************************************
!
!! SNRM2 computes the Euclidean norm of a vector.
!
!
!  Discussion:
!
!    The original SNRM2 algorithm is accurate but written in a bizarre,
!    unreadable and obsolete format.  This version goes for clarity.
!
!  Modified:
!
!    01 June 2000
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real X(*), the vector whose norm is to be computed.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real SNRM2, the Euclidean norm of X.
!
  integer i
  integer incx
  integer ix
  integer n
  real samax
  real snrm2
  real stemp
  real x(*)
  real xmax
!
  if ( n <= 0 ) then

    snrm2 = 0.0E+00

  else

    xmax = samax ( n, x, incx )

    if ( xmax == 0.0E+00 ) then

      snrm2 = 0.0E+00

    else

      if ( incx >= 0 ) then
        ix = 1
      else
        ix = ( - n + 1 ) * incx + 1
      end if

      stemp = 0.0E+00
      do i = 1, n
        stemp = stemp + ( x(ix) / xmax )**2
        ix = ix + incx
      end do

      snrm2 = xmax * sqrt ( stemp )

    end if

  end if

  return
end
subroutine snsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,ml, &
  mu,epsfcn,diag,mode,factor,nprint,info,nfev,njev,r,lr,qtf,wa1, &
  wa2,wa3,wa4)
!
!*******************************************************************************
!
!! SNSQ finds a zero of a system of n nonlinear functions in n variables.
!
!
!  by a modification of the powell
!            hybrid method.  this code is the combination of the minpack
!            codes (argonne) hybrd and hybrdj.
!
!  Description:
!
! 1. purpose.
!
!       the purpose of snsq is to find a zero of a system of n non-
!       linear functions in n variables by a modification of the powell
!       hybrid method.  the user must provide a subroutine which calcu-
!       lates the functions.  the user has the option of either to
!       provide a subroutine which calculates the jacobian or to let the
!       code calculate it by a forward-difference approximation.
!       this code is the combination of the minpack codes (argonne)
!       hybrd and hybrdj.
!
!
! 2. subroutine and type statements.
!
!     subroutine snsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,
!      *                 ml,mu,epsfcn,diag,mode,factor,nprint,info,nfev,
!      *                 njev,r,lr,qtf,wa1,wa2,wa3,wa4)
!       integer iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,njev,lr
!       real xtol,epsfcn,factor
!       real x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),qtf(n),
!      *     wa1(n),wa2(n),wa3(n),wa4(n)
!       external fcn,jac
!
!
! 3. parameters.
!
!       parameters designated as input parameters must be specified on
!       entry to snsq and are not changed on exit, while parameters
!       designated as output parameters need not be specified on entry
!       and are set to appropriate values on exit from snsq.
!
!       fcn is the name of the user-supplied subroutine which calculates
!         the functions.  fcn must be declared in an external statement
!         in the user calling program, and should be written as follows.
!
!       subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         real x(n),fvec(n)
!
!         calculate the functions at x and
!         return this vector in fvec.
!
!         return
!       end
!
!         the value of iflag should not be changed by fcn unless the
!         user wants to terminate execution of snsq.  in this case, set
!         iflag to a negative integer.
!
!       jac is the name of the user-supplied subroutine which calculates
!         the jacobian.  if iopt=1, then jac must be declared in an
!         external statement in the user calling program, and should be
!         written as follows.
!
!       subroutine jac(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         real x(n),fvec(n),fjac(ldfjac,n)
!
!         calculate the jacobian at x and return this
!         matrix in fjac.  fvec contains the function
!         values at x and should not be altered.
!
!         return
!       end
!
!         the value of iflag should not be changed by jac unless the
!         user wants to terminate execution of snsq.  in this case, set
!         iflag to a negative integer.
!
!         if iopt=2, jac can be ignored (treat it as a dummy argument).
!
!       iopt is an input variable which specifies how the jacobian will
!         be calculated.  if iopt=1, then the user must supply the
!         jacobian through the subroutine jac.  if iopt=2, then the
!         code will approximate the jacobian by forward-differencing.
!
!       n is a positive integer input variable set to the number of
!       functions and variables.
!
!       x is an array of length n.  on input, x must contain an initial
!         estimate of the solution vector.  on output, x contains the
!         final estimate of the solution vector.
!
!       fvec is an output array of length n which contains the functions
!         evaluated at the output x.
!
!       fjac is an output n by n array which contains the orthogonal
!         matrix q produced by the qr factorization of the final approx-
!         imate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       xtol is a non-negative input variable.  termination occurs when
!         the relative error between two consecutive iterates is at most
!         xtol.  therefore, xtol measures the relative error desired in
!         the approximate solution.  section 4 contains more details
!         about xtol.
!
!       maxfev is a positive integer input variable.  termination occurs
!         when the number of calls to fcn is at least maxfev by the end
!         of an iteration.
!
!       ml is a non-negative integer input variable which specifies the
!         number of subdiagonals within the band of the jacobian matrix.
!         if the jacobian is not banded or iopt=1, set ml to at
!         least n - 1.
!
!       mu is a non-negative integer input variable which specifies the
!         number of superdiagonals within the band of the jacobian
!         matrix.  if the jacobian is not banded or iopt=1, set mu to at
!         least n - 1.
!
!       epsfcn is an input variable used in determining a suitable step
!         for the forward-difference approximation.  this approximation
!         assumes that the relative errors in the functions are of the
!         order of epsfcn.  if epsfcn is less than the machine preci-
!         sion, it is assumed that the relative errors in the functions
!         are of the order of the machine precision.  if iopt=1, then
!         epsfcn can be ignored (treat it as a dummy argument).
!
!       diag is an array of length n.  if mode = 1 (see below), diag is
!         internally set.  if mode = 2, diag must contain positive
!         entries that serve as implicit (multiplicative) scale factors
!         for the variables.
!
!       mode is an integer input variable.  if mode = 1, the variables
!         will be scaled internally.  if mode = 2, the scaling is speci-
!         fied by the input diag.  other values of mode are equivalent
!         to mode = 1.
!
!       factor is a positive input variable used in determining the ini-
!         tial step bound.  this bound is set to the product of factor
!         and the euclidean norm of diag*x if nonzero, or else to factor
!         itself.  in most cases factor should lie in the interval
!         (.1,100.).  100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive.  in this case, fcn is
!         called with iflag = 0 at the beginning of the first iteration
!         and every nprint iteration thereafter and immediately prior
!         to return, with x and fvec available for printing. appropriate
!         print statements must be added to fcn(see example).  if nprint
!         is not positive, no special calls of fcn with iflag = 0 are
!         made.
!
!       info is an integer output variable.  if the user has terminated
!         execution, info is set to the (negative) value of iflag.  see
!         description of fcn and jac. otherwise, info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  relative error between two consecutive iterates is
!                   at most xtol.
!
!         info = 2  number of calls to fcn has reached or exceeded
!                   maxfev.
!
!         info = 3  xtol is too small.  no further improvement in the
!                   approximate solution x is possible.
!
!         info = 4  iteration is not making good progress, as measured
!                   by the improvement from the last five jacobian eval-
!                   uations.
!
!         info = 5  iteration is not making good progress, as measured
!                   by the improvement from the last ten iterations.
!
!         sections 4 and 5 contain more details about info.
!
!       nfev is an integer output variable set to the number of calls to
!         fcn.
!
!       njev is an integer output variable set to the number of calls to
!         jac. (if iopt=2, then njev is set to zero.)
!
!       r is an output array of length lr which contains the upper
!         triangular matrix produced by the qr factorization of the
!         final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains the vector
!         (q transpose)*fvec.
!
!       wa1, wa2, wa3, and wa4 are work arrays of length n.
!
!
! 4. successful completion.
!
!       the accuracy of snsq is controlled by the convergence parameter
!       xtol.  this parameter is used in a test which makes a comparison
!       between the approximation x and a solution xsol.  snsq termi-
!       nates when the test is satisfied.  if the convergence parameter
!       is less than the machine precision (as defined by the function
!       r1mach(4)), then snsq only attempts to satisfy the test
!       defined by the machine precision.  further progress is not
!       usually possible.
!
!       the test assumes that the functions are reasonably well behaved,
!       and, if the jacobian is supplied by the user, that the functions
!       and the jacobian are coded consistently.  if these conditions
!       are not satisfied, then snsq may incorrectly indicate conver-
!       gence.  the coding of the jacobian can be checked by the
!       subroutine chkder. if the jacobian is coded correctly or iopt=2,
!       then the validity of the answer can be checked, for example, by
!       rerunning snsq with a tighter tolerance.
!
!       convergence test.  if snrm2(z) denotes the euclidean norm of a
!         vector z and d is the diagonal matrix whose entries are
!         defined by the array diag, then this test attempts to guaran-
!         tee that
!
!               snrm2(d*(x-xsol)) <= xtol*snrm2(d*xsol).
!
!         if this condition is satisfied with xtol = 10**(-k), then the
!         larger components of d*x have k significant decimal digits and
!         info is set to 1.  there is a danger that the smaller compo-
!         nents of d*x may have large relative errors, but the fast rate
!         of convergence of snsq usually avoids this possibility.
!         unless high precision solutions are required, the recommended
!         value for xtol is the square root of the machine precision.
!
!
! 5. unsuccessful completion.
!
!       unsuccessful termination of snsq can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of func-
!       tion evaluations, or lack of good progress.
!
!       improper input parameters.  info is set to 0 if iopt < 1,
!         or iopt > 2, or n <= 0, or ldfjac < n, or
!         xtol < 0.0, or maxfev <= 0, or ml < 0, or mu < 0,
!         or factor <= 0.0, or lr < (n*(n+1))/2.
!
!       arithmetic interrupts.  if these interrupts occur in the fcn
!         subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of x by snsq.  in this
!         case, it may be possible to remedy the situation by rerunning
!         snsq with a smaller value of factor.
!
!       excessive number of function evaluations.  a reasonable value
!         for maxfev is 100*(n+1) for iopt=1 and 200*(n+1) for iopt=2.
!         if the number of calls to fcn reaches maxfev, then this
!         indicates that the routine is converging very slowly as
!         measured by the progress of fvec, and info is set to 2.  this
!         situation should be unusual because, as indicated below, lack
!         of good progress is usually diagnosed earlier by snsq,
!         causing termination with info = 4 or info = 5.
!
!       lack of good progress.  snsq searches for a zero of the system
!         by minimizing the sum of the squares of the functions.  in so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a zero of the system and, in this situ-
!         ation, the iteration eventually fails to make good progress.
!         in particular, this will happen if the system does not have a
!         zero.  if the system has a zero, rerunning snsq from a dif-
!         ferent starting point may be helpful.
!
!
! 6. characteristics of the algorithm.
!
!       snsq is a modification of the powell hybrid method.  two of its
!       main characteristics involve the choice of the correction as a
!       convex combination of the newton and scaled gradient directions,
!       and the updating of the jacobian by the rank-1 method of broy-
!       den.  the choice of the correction guarantees (under reasonable
!       conditions) global convergence for starting points far from the
!       solution and a fast rate of convergence.  the jacobian is
!       calculated at the starting point by either the user-supplied
!       subroutine or a forward-difference approximation, but it is not
!       recalculated until the rank-1 method fails to produce satis-
!       factory progress.
!
!       timing.  the time required by snsq to solve a given problem
!         depends on n, the behavior of the functions, the accuracy
!         requested, and the starting point.  the number of arithmetic
!         operations needed by snsq is about 11.5*(n**2) to process
!         each evaluation of the functions (call to fcn) and 1.3*(n**3)
!         to process each evaluation of the jacobian (call to jac,
!         if iopt = 1).  unless fcn and jac can be evaluated quickly,
!         the timing of snsq will be strongly influenced by the time
!         spent in fcn and jac.
!
!       storage.  snsq requires (3*n**2 + 17*n)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  there are no internally declared storage arrays.
!
!
! 7. example.
!
!       the problem is to determine the values of x(1), x(2), ..., x(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*x(1))*x(1)           -2*x(2)                   = -1
!               -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!                                   -x(8) + (3-2*x(9))*x(9) = -1
! c     **********
!
!       program test(input,output,tape6=output)
! c
! c     driver for snsq example.
! c
!       integer j,iopt,n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr,
!      *        nwrite
!       real xtol,epsfcn,factor,fnorm
!       real x(9),fvec(9),diag(9),fjac(9,9),r(45),qtf(9),
!      *     wa1(9),wa2(9),wa3(9),wa4(9)
!       real snrm2,r1mach
!       external fcn
!       data nwrite /6/
!
!       iopt = 2
!       n = 9
! c
! c     the following starting values provide a rough solution.
! c
!       do j = 1, 9
!          x(j) = -1.0E+00
!       end do
!
!       ldfjac = 9
!       lr = 45
! c
! c     set xtol to the square root of the machine precision.
! c     unless high precision solutions are required,
! c     this is the recommended setting.
! c
!       xtol = sqrt(r1mach(4))
! c
!       maxfev = 2000
!       ml = 1
!       mu = 1
!       epsfcn = 0.0E+00
!       mode = 2
!       do j = 1, 9
!          diag(j) = 1.0E+00
!       end do
!       factor = 1.e2
!       nprint = 0
! c
!       call snsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,ml,mu,
!      *           epsfcn,diag,mode,factor,nprint,info,nfev,njev,
!      *           r,lr,qtf,wa1,wa2,wa3,wa4)
!       fnorm = snrm2(n,fvec,1)
!       write (nwrite,1000) fnorm,nfev,info,(x(j),j=1,n)
!       stop
!  1000 format (5x,' final l2 norm of the residuals',e15.7 //
!      *        5x,' number of function evaluations',i10 //
!      *        5x,' exit parameter',16x,i10 //
!      *        5x,' final approximate solution' // (5x,3e15.7))
!       end
!       subroutine fcn(n,x,fvec,iflag)
!       integer n,iflag
!       real x(n),fvec(n)
!       integer k
!       real temp,temp1,temp2
! c
!       if (iflag /= 0) go to 5
! c
! c     insert print statements here when nprint is positive.
! c
!       return
!     5 continue
!       do k = 1, n
!          temp = (3.0E+00 - 2.0*x(k))*x(k)
!          temp1 = 0.0E+00
!          if (k /= 1) temp1 = x(k-1)
!          temp2 = 0.0E+00
!          if (k /= n) temp2 = x(k+1)
!          fvec(k) = temp - temp1 - 2.0*temp2 + 1.0E+00
!       end do
!       return
!       end
!
!       results obtained with different compilers or machines
!       may be slightly different.
!
!       final l2 norm of the residuals  0.1192636e-07
!
!       number of function evaluations        14
!
!       exit parameter                         1
!
!       final approximate solution
!
!       -0.5706545e+00 -0.6816283e+00 -0.7017325e+00
!       -0.7042129e+00 -0.7013690e+00 -0.6918656e+00
!       -0.6657920e+00 -0.5960342e+00 -0.4164121e+00
!
!  references  powell, m. j. d.
!                 a hybrid method for nonlinear equations.
!                 numerical methods for nonlinear algebraic equations,
!                 p. rabinowitz, editor.  gordon and breach, 1970.
!
  integer ldfjac
  integer lr
  integer n
!
  real actred
  real delta
  real diag(n)
  real epsfcn
  real epsmch
  real factor
  real fjac(ldfjac,n)
  real fnorm
  real fnorm1
  real fvec(n)
  integer i
  integer iflag
  integer info
  integer iopt
  integer iter
  integer iwa(1)
  integer j
  logical jeval
  integer jm1
  integer l
  integer maxfev
  integer ml
  integer mode
  integer mu
  integer ncfail
  integer ncsuc
  integer nfev
  integer njev
  integer nprint
  integer nslow1
  integer nslow2
  real p001
  real p0001
  real, parameter :: p1 = 0.1E+00
  real, parameter :: p5 = 0.5E+00
  real pnorm
  real prered
  real qtf(n)
  real r(lr)
  real r1mach
  real ratio
  logical sing
  real snrm2
  real sum
  real temp
  real xnorm
  real wa1(n)
  real wa2(n)
  real wa3(n)
  real wa4(n)
  real x(n)
  real xtol
!
  external fcn
  external jac
!
  data p001,p0001 /1.0e-03,1.0e-04/
!
  epsmch = epsilon ( epsmch )
  info = 0
  iflag = 0
  nfev = 0
  njev = 0

  if (iopt < 1 .or. iopt > 2 ) then
    go to 300
  else if ( n <= 0 .or. xtol < 0.0E+00 .or. maxfev <= 0 ) then
    go to 300
  else if ( ml < 0 .or. mu < 0 .or. factor <= 0.0E+00 ) then
    go to 300
  else if ( ldfjac < n .or. lr < (n*(n + 1))/2) then
    go to 300
  end if

  if ( mode /= 2 ) then
    go to 20
  end if

  do j = 1, n
    if ( diag(j) <= 0.0E+00 ) then
      go to 300
    end if
  end do

   20 continue
!
!  evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn ( n, x, fvec, iflag )
  nfev = 1
  if (iflag < 0) go to 300
  fnorm = snrm2 ( n, fvec, 1 )
!
!  initialize iteration counter and monitors.
!
  iter = 1
  ncsuc = 0
  ncfail = 0
  nslow1 = 0
  nslow2 = 0
!
!  beginning of the outer loop.
!
   30 continue
     jeval = .true.
!
!  calculate the jacobian matrix.
!
     if ( iopt /= 2 ) then
!
!  user supplies jacobian
!
        call jac(n,x,fvec,fjac,ldfjac,iflag)
        njev = njev+1
!
!  code approximates the jacobian
!
      else

        iflag = 2

        call fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )

        nfev = nfev + min(ml+mu+1,n)

     end if

     if (iflag < 0) go to 300
!
!  compute the qr factorization of the jacobian.
!
     call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!  on the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if (iter /= 1) go to 70
     if (mode == 2) go to 50
     diag(1:n) = wa2(1:n)
     do j = 1, n
       if (wa2(j) == 0.0) diag(j) = 1.0E+00
     end do

   50    continue
!
!  on the first iteration, calculate the norm of the scaled x
!  and initialize the step bound delta.
!
     wa3(1:n) = diag(1:n) * x(1:n)
     xnorm = snrm2(n,wa3,1)
     delta = factor*xnorm
     if (delta == 0.0) delta = factor

70   continue
!
!  form (q transpose)*fvec and store in qtf.
!
     qtf(1:n) = fvec(1:n)

     do j = 1, n

       if ( fjac(j,j) /= 0.0E+00 ) then

        sum = 0.0E+00
        do i = j, n
           sum = sum + fjac(i,j)*qtf(i)
        end do

        temp = -sum/fjac(j,j)
        do i = j, n
           qtf(i) = qtf(i) + fjac(i,j)*temp
        end do

      end if

    end do
!
!  copy the triangular factor of the qr factorization into r.
!
     sing = .false.

     do j = 1, n

        l = j
        jm1 = j - 1
        do i = 1, jm1
           r(l) = fjac(i,j)
           l = l + n - i
        end do
        r(l) = wa1(j)
        if (wa1(j) == 0.0) sing = .true.

     end do
!
!  accumulate the orthogonal factor in fjac.
!
     call qform(n,n,fjac,ldfjac,wa1)
!
!  rescale if necessary.
!
     if (mode == 2) go to 170

     do j = 1, n
       diag(j) = max (diag(j),wa2(j))
     end do

  170    continue
!
!  beginning of the inner loop.
!
  180    continue

!
!  if requested, call fcn to enable printing of iterates.
!
        if (nprint > 0) then
          iflag = 0
          if (mod(iter-1,nprint) == 0) call fcn(n,x,fvec,iflag)
          if (iflag < 0) go to 300
        end if
!
!  determine the direction p.
!
        call dogleg ( n, r, lr, diag, qtf, delta, wa1, wa2, wa3 )

!
!  store the direction p and x + p. calculate the norm of p.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n)*wa1(1:n)
        pnorm = snrm2(n,wa3,1)
!
!  on the first iteration, adjust the initial step bound.
!
        if (iter == 1) delta = min (delta,pnorm)
!
!  Evaluate the function at x + p and calculate its norm.
!
        iflag = 1
        call fcn(n,wa2,wa4,iflag)
        nfev = nfev + 1
        if (iflag < 0) go to 300
        fnorm1 = snrm2(n,wa4,1)
!
!  Compute the scaled actual reduction.
!
        actred = - 1.0E+00
        if (fnorm1 < fnorm) actred = 1.0E+00 - (fnorm1/fnorm)**2
!
!  compute the scaled predicted reduction.
!
        l = 1
        do i = 1, n
           sum = 0.0E+00
           do j = i, n
              sum = sum + r(l)*wa1(j)
              l = l + 1
           end do
           wa3(i) = qtf(i) + sum
        end do

        temp = snrm2(n,wa3,1)
        prered = 0.0E+00
        if (temp < fnorm) prered = 1.0E+00 - (temp/fnorm)**2
!
!  Compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0E+00
        if (prered > 0.0) ratio = actred/prered
!
!  update the step bound.
!
        if ( ratio < p1 ) then

           ncsuc = 0
           ncfail = ncfail + 1
           delta = p5*delta

        else

           ncfail = 0
           ncsuc = ncsuc + 1

           if (ratio >= p5 .or. ncsuc > 1) then
             delta = max (delta,pnorm/p5)
           end if

           if (abs(ratio-1.0) <= p1) delta = pnorm/p5

        end if
!
!  successful iteration. update x, fvec, and their norms.
!
        if ( ratio >= p0001 ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n)*x(1:n)
          fvec(1:n) = wa4(1:n)
          xnorm = snrm2(n,wa2,1)
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  determine the progress of the iteration.
!
        nslow1 = nslow1 + 1
        if (actred >= p001) nslow1 = 0
        if (jeval) nslow2 = nslow2 + 1
        if (actred >= p1) nslow2 = 0
!
!  test for convergence.
!
        if (delta <= xtol*xnorm .or. fnorm == 0.0) info = 1
        if (info /= 0) go to 300
!
!  tests for termination and stringent tolerances.
!
        if (nfev >= maxfev) info = 2
        if (p1* max (p1*delta,pnorm) <= epsmch*xnorm) info = 3
        if (nslow2 == 5) info = 4
        if (nslow1 == 10) info = 5
        if (info /= 0) go to 300
!
!  criterion for recalculating jacobian
!
        if (ncfail == 2) go to 290
!
!  calculate the rank one modification to the jacobian
!  and update qtf if necessary.
!
        do j = 1, n
           sum = 0.0E+00
           do i = 1, n
              sum = sum + fjac(i,j)*wa4(i)
           end do
           wa2(j) = (sum - wa3(j))/pnorm
           wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
           if (ratio >= p0001) qtf(j) = sum
        end do
!
!  compute the qr factorization of the updated jacobian.
!
        call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
        call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
        call r1mpyq(1,n,qtf,1,wa2,wa3)
!
!  end of the inner loop.
!
        jeval = .false.
        go to 180
  290    continue
!
!  end of the outer loop.
!
     go to 30
  300 continue
!
!  termination, either normal or user imposed.
!
  if (iflag < 0) info = iflag
  iflag = 0
  if (nprint > 0) call fcn(n,x,fvec,iflag)
  if (info < 0) call xerror( &
    'snsq   -- execution terminated because user set iflag negative.',63,1,1)

  if (info == 0) call xerror( 'snsq   -- invalid input parameter.',34,2,1)
  if (info == 2) call xerror( 'snsq   -- too many function evaluations.',40,9,1)
  if (info == 3) then
    call xerror( 'snsq   -- xtol too small. no further improvement possible.', &
      58,3,1)
  end if

  if ( info > 4 ) then
    call xerror( 'snsq   -- iteration not making good progress.',45,1,1)
  end if

  return
end
subroutine snsqe(fcn,jac,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
!
!*******************************************************************************
!
!! SNSQE is the easy-to-use version of SNSQ.
!
!
!  Discussion:
!
!    SNSQ finds a zero of a system of N nonlinear functions in N variables by a
!    modification of the Powell hybrid method.  This code is the
!    combination of the MINPACK codes HYBRD1 and HYBRJ1.
!
!  Description:
!
! 1. purpose.
!
!
!       the purpose of snsqe is to find a zero of a system of n non-
!       linear functions in n variables by a modification of the powell
!       hybrid method.  this is done by using the more general nonlinear
!       equation solver snsq.  the user must provide a subroutine which
!       calculates the functions.  the user has the option of either to
!       provide a subroutine which calculates the jacobian or to let the
!       code calculate it by a forward-difference approximation.  this
!       code is the combination of the minpack codes (argonne) hybrd1
!       and hybrj1.
!
!
! 2. subroutine and type statements.
!
!     subroutine snsqe(fcn,jac,iopt,n,x,fvec,tol,nprint,info,
!      *                  wa,lwa)
!       integer iopt,n,nprint,info,lwa
!       real tol
!       real x(n),fvec(n),wa(lwa)
!       external fcn,jac
!
! 3. parameters.
!
!       parameters designated as input parameters must be specified on
!       entry to snsqe and are not changed on exit, while parameters
!       designated as output parameters need not be specified on entry
!       and are set to appropriate values on exit from snsqe.
!
!       fcn is the name of the user-supplied subroutine which calculates
!         the functions.  fcn must be declared in an external statement
!         in the user calling program, and should be written as follows.
!
!       subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         real x(n),fvec(n)
!
!         calculate the functions at x and
!         return this vector in fvec.
!
!         return
!       end
!
!         the value of iflag should not be changed by fcn unless the
!         user wants to terminate execution of snsqe.  in this case, set
!         iflag to a negative integer.
!
!       jac is the name of the user-supplied subroutine which calculates
!         the jacobian.  if iopt=1, then jac must be declared in an
!         external statement in the user calling program, and should be
!         written as follows.
!
!       subroutine jac(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         real x(n),fvec(n),fjac(ldfjac,n)
!
!         calculate the jacobian at x and return this
!         matrix in fjac.  fvec contains the function
!         values at x and should not be altered.
!
!         return
!       end
!
!         the value of iflag should not be changed by jac unless the
!         user wants to terminate execution of snsqe.  in this case, set
!         iflag to a negative integer.
!
!         if iopt=2, jac can be ignored (treat it as a dummy argument).
!
!       iopt is an input variable which specifies how the jacobian will
!         be calculated.  if iopt=1, then the user must supply the
!         jacobian through the subroutine jac.  if iopt=2, then the
!         code will approximate the jacobian by forward-differencing.
!
!       n is a positive integer input variable set to the number of
!       functions and variables.
!
!       x is an array of length n.  on input, x must contain an initial
!         estimate of the solution vector.  on output, x contains the
!         final estimate of the solution vector.
!
!       fvec is an output array of length n which contains the functions
!         evaluated at the output x.
!
!       tol is a non-negative input variable.  termination occurs when
!         the algorithm estimates that the relative error between x and
!         the solution is at most tol.  section 4 contains more details
!         about tol.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive.  in this case, fcn is
!         called with iflag = 0 at the beginning of the first iteration
!         and every nprint iteration thereafter and immediately prior
!         to return, with x and fvec available for printing. appropriate
!         print statements must be added to fcn (see example). if nprint
!         is not positive, no special calls of fcn with iflag = 0 are
!         made.
!
!       info is an integer output variable.  if the user has terminated
!         execution, info is set to the (negative) value of iflag.  see
!         description of fcn and jac. otherwise, info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error between
!                   x and the solution is at most tol.
!
!         info = 2  number of calls to fcn has reached or exceeded
!                   100*(n+1) for iopt=1 or 200*(n+1) for iopt=2.
!
!         info = 3  tol is too small.  no further improvement in the
!                   approximate solution x is possible.
!
!         info = 4  iteration is not making good progress.
!
!         sections 4 and 5 contain more details about info.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         (3*n**2+13*n))/2.
!
!
! 4. successful completion.
!
!       the accuracy of snsqe is controlled by the convergence parame-
!       ter tol.  this parameter is used in a test which makes a compar-
!       ison between the approximation x and a solution xsol.  snsqe
!       terminates when the test is satisfied.  if tol is less than the
!       machine precision (as defined by the function r1mach(4)), then
!       snsqe attemps only to satisfy the test defined by the machine
!       precision.  further progress is not usually possible.  unless
!       high precision solutions are required, the recommended value
!       for tol is the square root of the machine precision.
!
!       the test assumes that the functions are reasonably well behaved,
!       and, if the jacobian is supplied by the user, that the functions
!       and the jacobian  coded consistently.  if these conditions
!       are not satisfied, snsqe may incorrectly indicate convergence.
!       the coding of the jacobian can be checked by the subroutine
!       chkder.  if the jacobian is coded correctly or iopt=2, then
!       the validity of the answer can be checked, for example, by
!       rerunning snsqe with a tighter tolerance.
!
!       convergence test.  if snrm2(z) denotes the euclidean norm of a
!         vector z, then this test attempts to guarantee that
!
!               snrm2(x-xsol) <=  tol*snrm2(xsol).
!
!         if this condition is satisfied with tol = 10**(-k), then the
!         larger components of x have k significant decimal digits and
!         info is set to 1.  there is a danger that the smaller compo-
!         nents of x may have large relative errors, but the fast rate
!         of convergence of snsqe usually avoids this possibility.
!
!
! 5. unsuccessful completion.
!
!       unsuccessful termination of snsqe can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of func-
!       tion evaluations, errors in the functions, or lack of good prog-
!       ress.
!
!       improper input parameters.  info is set to 0 if iopt < 1, or
!         iopt > 2, or n <= 0, or tol < 0.0, or
!         lwa < (3*n**2+13*n)/2.
!
!       arithmetic interrupts.  if these interrupts occur in the fcn
!       subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of x by snsqe.  in this
!         case, it may be possible to remedy the situation by not evalu-
!         ating the functions here, but instead setting the components
!         of fvec to numbers that exceed those in the initial fvec.
!
!       excessive number of function evaluations.  if the number of
!         calls to fcn reaches 100*(n+1) for iopt=1 or 200*(n+1) for
!         iopt=2, then this indicates that the routine is converging
!         very slowly as measured by the progress of fvec, and info is
!         set to 2.  this situation should be unusual because, as
!         indicated below, lack of good progress is usually diagnosed
!         earlier by snsqe, causing termination with info = 4.
!
!       errors in the functions.  when iopt=2, the choice of step length
!         in the forward-difference approximation to the jacobian
!         assumes that the relative errors in the functions are of the
!         order of the machine precision.  if this is not the case,
!         snsqe may fail (usually with info = 4).  the user should
!         then either use snsq and set the step length or use iopt=1
!         and supply the jacobian.
!
!       lack of good progress.  snsqe searches for a zero of the system
!         by minimizing the sum of the squares of the functions.  in so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a zero of the system and, in this situ-
!         ation, the iteration eventually fails to make good progress.
!         in particular, this will happen if the system does not have a
!         zero.  if the system has a zero, rerunning snsqe from a dif-
!         ferent starting point may be helpful.
!
!
! 6. characteristics of the algorithm.
!
!       snsqe is a modification of the powell hybrid method.  two of
!       its main characteristics involve the choice of the correction as
!       a convex combination of the newton and scaled gradient direc-
!       tions, and the updating of the jacobian by the rank-1 method of
!       broyden.  the choice of the correction guarantees (under reason-
!       able conditions) global convergence for starting points far from
!       the solution and a fast rate of convergence.  the jacobian is
!       calculated at the starting point by either the user-supplied
!       subroutine or a forward-difference approximation, but it is not
!       recalculated until the rank-1 method fails to produce satis-
!       factory progress.
!
!       timing.  the time required by snsqe to solve a given problem
!         depends on n, the behavior of the functions, the accuracy
!         requested, and the starting point.  the number of arithmetic
!         operations needed by snsqe is about 11.5*(n**2) to process
!         each evaluation of the functions (call to fcn) and 1.3*(n**3)
!         to process each evaluation of the jacobian (call to jac,
!         if iopt = 1).  unless fcn and jac can be evaluated quickly,
!         the timing of snsqe will be strongly influenced by the time
!         spent in fcn and jac.
!
!       storage.  snsqe requires (3*n**2 + 17*n)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  there are no internally declared storage arrays.
!
!
! 7. example.
!
!       the problem is to determine the values of x(1), x(2), ..., x(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*x(1))*x(1)           -2*x(2)                   = -1
!               -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!                                   -x(8) + (3-2*x(9))*x(9) = -1
!
!       program test(input,output,tape6=output)
! c
! c     driver for snsqe example.
! c
!       integer j,n,iopt,nprint,info,lwa,nwrite
!       real tol,fnorm
!       real x(9),fvec(9),wa(180)
!       real snrm2,r1mach
!       external fcn
!       data nwrite /6/
! c
!       iopt = 2
!       n = 9
! c
! c     the following starting values provide a rough solution.
! c
!       do j = 1, 9
!          x(j) = -1.0E+00
!       end do
!
!       lwa = 180
!       nprint = 0
! c
! c     set tol to the square root of the machine precision.
! c     unless high precision solutions are required,
! c     this is the recommended setting.
! c
!       tol = sqrt(r1mach(4))
! c
!       call snsqe(fcn,jac,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
!       fnorm = snrm2(n,fvec)
!       write (nwrite,1000) fnorm,info,(x(j),j=1,n)
!       stop
!  1000 format (5x,' final l2 norm of the residuals',e15.7 //
!      *        5x,' exit parameter',16x,i10 //
!      *        5x,' final approximate solution' // (5x,3e15.7))
!     end
!     subroutine fcn(n,x,fvec,iflag)
!       integer n,iflag
!       real x(n),fvec(n)
!       integer k
!       real temp,temp1,temp2
!
!       do k = 1, n
!          temp = (3.0E+00 - 2.0*x(k))*x(k)
!          temp1 = 0.0E+00
!          if (k /= 1) temp1 = x(k-1)
!          temp2 = 0.0E+00
!          if (k /= n) temp2 = x(k+1)
!          fvec(k) = temp - temp1 - 2.0*temp2 + 1.0E+00
!       end do
!
!       return
!     end
!
!       results obtained with different compilers or machines
!       may be slightly different.
!
!       final l2 norm of the residuals  0.1192636e-07
!
!       exit parameter                         1
!
!       final approximate solution
!
!       -0.5706545e+00 -0.6816283e+00 -0.7017325e+00
!       -0.7042129e+00 -0.7013690e+00 -0.6918656e+00
!       -0.6657920e+00 -0.5960342e+00 -0.4164121e+00
!
!  references  powell, m. j. d.
!                 a hybrid method for nonlinear equations.
!                 numerical methods for nonlinear algebraic equations,
!                 p. rabinowitz, editor.  gordon and breach, 1970.
!
  integer lwa
  integer n
!
  real epsfcn
  real, parameter :: factor = 100.0E+00
  real fvec(n)
  integer index
  integer info
  integer iopt
  integer j
  integer lr
  integer maxfev
  integer ml
  integer mode
  integer mu
  integer nfev
  integer njev
  integer nprint
  real tol
  real wa(lwa)
  real x(n)
  real xtol
!
  external fcn
  external jac
!
  info = 0
!
!  Check the input parameters for errors.
!
  if (iopt < 1 .or. iopt > 2 .or. n <= 0 .or. tol < 0.0E+00 .or. &
    lwa < (3*n**2 +13*n)/2) then
    go to 20
  end if

  maxfev = 100*(n + 1)
  if ( iopt == 2 ) maxfev = 2*maxfev
  xtol = tol
  ml = n - 1
  mu = n - 1
  epsfcn = 0.0E+00
  mode = 2
  wa(1:n) = 1.0E+00
  lr = (n*(n + 1))/2
  index = 6*n + lr

  call snsq ( fcn, jac, iopt, n, x, fvec, wa(index+1), n, xtol, maxfev, ml, mu, &
    epsfcn, wa(1), mode, factor, nprint, info, nfev, njev, &
    wa(6*n+1), lr, wa(n+1), wa(2*n+1), wa(3*n+1), wa(4*n+1), wa(5*n+1) )

  if (info == 5) info = 4
   20 continue

  if (info == 0) then
    call xerror( 'snsqe  -- invalid input parameter.',34,2,1)
  end if

  return
end
subroutine sqrank(a,lda,m,n,tol,kr,jpvt,qraux,work)
!
!*******************************************************************************
!
!! SQRANK computes the QR factorization of a rectangular matrix.
!
!
!    this routine is used in conjunction with sqrlss to solve linear
!    systems of equations in a least square sense.
!
!  Description:
!
!     sqrank is used in conjunction with sqrlss to solve
!     overdetermined, underdetermined and singular linear systems
!     in a least squares sense.
!     sqrank uses the linpack subroutine sqrdc to compute the qr
!     factorization, with column pivoting, of an  m  by  n  matrix  a .
!     the numerical rank is determined using the tolerance tol.
!
!     for more information, see chapter 9 of linpack users guide,
!     j. dongarra et all, siam, 1979.
!
!     on entry
!
!        a     real (lda,n) .
!              the matrix whose decomposition is to be computed.
!
!        lda   integer.
!              the leading dimension of a .
!
!        m     integer.
!              the number of rows of a .
!
!        n     integer.
!              the number of columns of  a .
!
!        tol   real.
!              a relative tolerance used to determine the numerical
!              rank.  the problem should be scaled so that all the elements
!              of  a   have roughly the same absolute accuracy, eps.  then a
!              reasonable value for  tol  is roughly  eps  divided by
!              the magnitude of the largest element.
!
!        jpvt  integer(n)
!
!        qraux real(n)
!
!        work  real(n)
!              three auxilliary vectors used by sqrdc .
!
!     on return
!
!        a     contains the output from sqrdc.
!              the triangular matrix  r  of the qr factorization is
!              contained in the upper triangle and information needed
!              to recover the orthogonal matrix  q  is stored below
!              the diagonal in  a  and in the vector qraux .
!
!        kr    integer.
!              the numerical rank.
!
!        jpvt  the pivot information from sqrdc.
!
!     columns jpvt(1),...,jpvt(kr) of the original matrix are linearly
!     independent to within the tolerance tol and the remaining columns
!     are linearly dependent.  abs(a(1,1))/abs(a(kr,kr))  is an estimate
!     of the condition number of the matrix of independent columns,
!     and of  r .  this estimate will be <= 1/tol .
!
!      usage.....see subroutine sqrlss
!
!  reference(s)
!
!      dongarra, et al, linpack users guide, siam, 1979
!
  integer lda
  integer n
!
  real a(lda,n)
  integer j
  integer jpvt(*)
  integer k
  integer kr
  integer m
  real qraux(*)
  real tol
  real work(*)
!
  jpvt(1:n) = 0

  call sqrdc(a,lda,m,n,qraux,jpvt,work,1)
  kr = 0
  k = min(m,n)

  do j = 1, k
    if (abs(a(j,j)) <= tol*abs(a(1,1))) then
      return
    end if
    kr = j
  end do

  return
end
subroutine sqrdc(x,ldx,n,p,qraux,jpvt,work,job)
!
!*******************************************************************************
!
!! SQRDC computes the QR factorization of a rectangular matrix.
!
!
!  column pivoting
!     based on the 2-norms of the reduced columns may be
!     performed at the user's option.
!
!     on entry
!
!        x       real(ldx,p), where ldx >= n.
!                x contains the matrix whose decomposition is to be
!                computed.
!
!        ldx     integer.
!                ldx is the leading dimension of the array x.
!
!        n       integer.
!                n is the number of rows of the matrix x.
!
!        p       integer.
!                p is the number of columns of the matrix x.
!
!        jpvt    integer(p).
!                jpvt contains integers that control the selection
!                of the pivot columns.  the k-th column x(k) of x
!                is placed in one of three classes according to the
!                value of jpvt(k).
!
!                   if jpvt(k) > 0, then x(k) is an initial
!                                      column.
!
!                   if jpvt(k) == 0, then x(k) is a free column.
!
!                   if jpvt(k) < 0, then x(k) is a final column.
!
!                before the decomposition is computed, initial columns
!                are moved to the beginning of the array x and final
!                columns to the end.  both initial and final columns
!                are frozen in place during the computation and only
!                free columns are moved.  at the k-th stage of the
!                reduction, if x(k) is occupied by a free column,
!                it is interchanged with the free column of largest
!                reduced norm.  jpvt is not referenced if
!                job == 0.
!
!        work    real(p).
!                work is a work array.  work is not referenced if
!                job == 0.
!
!        job     integer.
!                job is an integer that initiates column pivoting.
!                if job == 0, no pivoting is done.
!                if job /= 0, pivoting is done.
!
!     on return
!
!        x       x contains in its upper triangle the upper
!                triangular matrix r of the qr factorization.
!                below its diagonal x contains information from
!                which the orthogonal part of the decomposition
!                can be recovered.  note that if pivoting has
!                been requested, the decomposition is not that
!                of the original matrix x but that of x
!                with its columns permuted as described by jpvt.
!
!        qraux   real(p).
!                qraux contains further information required to recover
!                the orthogonal part of the decomposition.
!
!        jpvt    jpvt(k) contains the index of the column of the
!                original matrix that has been interchanged into
!                the k-th column, if pivoting was requested.
!
!     linpack.  this version dated 08/14/78 .
!     g. w. stewart, university of maryland, argonne national lab.
!
!  references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
!                 *linpack users  guide*, siam, 1979.
!
  integer ldx
  integer p
!
  integer j
  integer jj
  integer job
  integer jp
  integer jpvt(p)
  integer l
  integer lp1
  integer lup
  integer maxj
  real maxnrm
  integer n
  logical negj
  real nrmxl
  integer pl
  integer pu
  real qraux(1)
  real sdot
  real snrm2
  logical swapj
  real t
  real tt
  real work(1)
  real x(ldx,1)
!
  pl = 1
  pu = 0

  if ( job /= 0 ) then
!
!  pivoting has been requested.  rearrange the columns according to jpvt.
!
     do j = 1, p
        swapj = jpvt(j) > 0
        negj = jpvt(j) < 0
        jpvt(j) = j
        if (negj) jpvt(j) = -j
        if (swapj) then
           if (j /= pl) call sswap(n,x(1,pl),1,x(1,j),1)
           jpvt(j) = jpvt(pl)
           jpvt(pl) = j
           pl = pl + 1
        end if
     end do

     pu = p

     do jj = 1, p
        j = p - jj + 1
        if (jpvt(j) < 0) then
           jpvt(j) = -jpvt(j)
           if (j /= pu) then
              call sswap(n,x(1,pu),1,x(1,j),1)
              jp = jpvt(pu)
              jpvt(pu) = jpvt(j)
              jpvt(j) = jp
           end if
           pu = pu - 1
        end if
      end do

  end if
!
!  compute the norms of the free columns.
!
  do j = pl, pu
     qraux(j) = snrm2(n,x(1,j),1)
     work(j) = qraux(j)
  end do
!
!  perform the householder reduction of x.
!
  lup = min(n,p)

  do l = 1, lup

     if ( pl <= l .and. l < pu ) then
!
!  locate the column of largest norm and bring it into the pivot position.
!
        maxnrm = 0.0E+00
        maxj = l

        do j = l, pu
          if ( qraux(j) > maxnrm ) then
            maxnrm = qraux(j)
            maxj = j
          end if
        end do

        if (maxj /= l) then
           call sswap(n,x(1,l),1,x(1,maxj),1)
           qraux(maxj) = qraux(l)
           work(maxj) = work(l)
           jp = jpvt(maxj)
           jpvt(maxj) = jpvt(l)
           jpvt(l) = jp
        end if

     end if

     qraux(l) = 0.0E+00

     if (l == n) go to 190
!
!  compute the householder transformation for column l.
!
        nrmxl = snrm2(n-l+1,x(l,l),1)

        if (nrmxl == 0.0) go to 180

           if (x(l,l) /= 0.0) nrmxl = sign(nrmxl,x(l,l))
           call sscal(n-l+1,1.0/nrmxl,x(l,l),1)
           x(l,l) = 1.0E+00 + x(l,l)
!
!  apply the transformation to the remaining columns, updating the norms.
!
           lp1 = l + 1

           do j = lp1, p

              t = -sdot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
              call saxpy(n-l+1,t,x(l,l),1,x(l,j),1)

              if (j < pl .or. j > pu) go to 150

              if ( qraux(j) /= 0.0E+00 ) then

                 tt = 1.0E+00 - (abs(x(l,j))/qraux(j))**2
                 tt = max (tt,0.0)
                 t = tt
                 tt = 1.0E+00 + 0.05*tt*(qraux(j)/work(j))**2

                 if (tt /= 1.0) then
                    qraux(j) = qraux(j)*sqrt(t)
                 else
                    qraux(j) = snrm2(n-l,x(l+1,j),1)
                    work(j) = qraux(j)
                 end if

               end if

  150             continue
           end do
!
!  save the transformation.
!
           qraux(l) = x(l,l)
           x(l,l) = -nrmxl
  180       continue
  190    continue

  end do

  return
end
subroutine sqrls (a,lda,m,n,tol,kr,b,x,rsd,jpvt,qraux,work,itask,ind)
!
!*******************************************************************************
!
!! SQRLS solves an linear system in the least squares sense.
!
!
!  Discussion:
!
!    The linear system may be overdetermined, underdetermined or singular.
!    The solution is obtained using a QR factorization of the
!    coefficient matrix.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     sqrls is used to solve in a least squares sense
!     overdetermined, underdetermined and singular linear systems .
!     the system is  a*x  approximates  b  where  a  is  m  by  n.
!     b  is a given  m-vector, and  x  is the  n-vector to be computed.
!     a solution  x  is found which minimimzes the sum of squares (2-norm)
!     of the residual,  a*x - b .
!     the numerical rank of a is determined using the tolerance tol.
!
!     sqrls uses the linpack subroutine sqrdc to compute the qr
!     factorization, with column pivoting, of an  m  by  n  matrix  a .
!     for more information, see chapter 9 of the reference below.
!
!
!     on entry
!
!        a     real (lda,n) .
!              the matrix whose decomposition is to be computed.
!              in a least squares data fitting problem, a(i,j) is the
!              value of the j-th basis (model) function at the i-th
!              data point.
!
!        lda   integer.
!              the leading dimension of a .
!
!        m     integer.
!              the number of rows of a .
!
!        n     integer.
!              the number of columns of  a .
!
!        tol   real.
!              a relative tolerance used to determine the numerical
!              rank.  the problem should be scaled so that all the
!              elements of  a   have roughly the same absolute accuracy
!              eps.  then a reasonable value for  tol  is roughly  eps
!              divided by the magnitude of the largest element.
!
!        jpvt  integer(n)
!        qraux real(n)
!        work  real(n)
!                    three auxiliary arrays used to factor the matrix a.
!                    (not required if itask > 1)
!
!        b     real(m)
!              the right hand side of the linear system.
!              in a least squares data fitting problem b(i) contains the
!              value of i-th observation.
!
!        itask integer.
!              if itask=1, then sqrls factors the matrix a and
!                          solves the least squares problem.
!              if itask=2, then sqrls assumes that the matrix a
!                          was factored with an earlier call to
!                          sqrls, and only solves the least squares
!                          problem.
!
!     on return
!
!        x     real(n) .
!              a least squares solution to the linear system.
!
!        rsd   real(m) .
!              the residual, b - a*x .  rsd may overwrite  b .
!
!        ind   integer
!              error code:  ind = 0:  no error
!                           ind = -1: n > lda   (fatal error)
!                           ind = -2: n < 1     (fatal error)
!                           ind = -3: itask < 1 (fatal error)
!
!        a     contains the output from sqrdc.
!              the triangular matrix  r  of the qr factorization is
!              contained in the upper triangle and information needed
!              to recover the orthogonal matrix  q  is stored below
!              the diagonal in  a  and in the vector qraux .
!
!        kr    integer.
!              the numerical rank.
!
!        jpvt  the pivot information from sqrdc.
!
!     columns jpvt(1),...,jpvt(kr) of the original matrix are linearly
!     independent to within the tolerance tol and the remaining columns
!     are linearly dependent.  abs(a(1,1))/abs(a(kr,kr))  is an estimate
!     of the condition number of the matrix of independent columns,
!     and of  r .  this estimate will be <= 1/tol .
!
!      usage....
!        sqrls can be efficiently used to solve several least squares
!      problems with the same matrix a.  the first system is solved
!      with itask = 1.  the subsequent systems are solved with
!      itask = 2, to avoid the recomputation of the matrix factors.
!      the parameters  kr, jpvt, and qraux must not be modified
!      between calls to sqrls.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real b(m)
  integer ind
  integer itask
  integer jpvt(n)
  integer kr
  real qraux(n)
  real rsd(m)
  real tol
  real work(n)
  real x(n)
!
  if ( lda < n ) then
    write ( *, * ) ' '
    write ( *, * ) 'SQRLS - Fatal error!'
    write ( *, * ) '  LDA < N.'
    stop
  end if

  if ( n <= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SQRLS - Fatal error!'
    write ( *, * ) '  N <= 0.'
    stop
  end if

  if ( itask < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SQRLS - Fatal error!'
    write ( *, * ) '  ITASK < 1.'
    stop
  end if

  ind = 0
!
!  Factor the matrix.
!
  if ( itask == 1 ) then
    call sqrank ( a, lda, m, n, tol, kr, jpvt, qraux, work )
  end if
!
!  Solve the least-squares problem.
!
  call sqrlss ( a, lda, m, n, kr, b, x, rsd, jpvt, qraux )

  return
end
subroutine sqrlss ( a, lda, m, n, kr, b, x, rsd, jpvt, qraux )
!
!*******************************************************************************
!
!! SQRLSS solves a linear system in a least squares sense.
!
!
!  Discussion:
!
!    SQRLSS must be preceeded by a call to SQRANK.
!
!    The system is to be solved is
!      A * X = B
!    where
!      A is an M by N matrix with rank KR, as determined by SQRANK,
!      B is a given M-vector,
!      X is the N-vector to be computed.
!
!    A solution X, with at most KR nonzero components, is found which
!    minimizes the 2-norm of the residual (A*X-B).
!
!     on entry
!
!        a,lda,m,n,kr,jpvt,qraux
!              the output from sqrank .
!
!        b     real(m) .
!              the right hand side of the linear system.
!
!     on return
!
!        x     real(n) .
!              a least squares solution to the linear system.
!
!        rsd   real(m) .
!              the residual, b - a*x .  rsd may overwite  b .
!
!      usage....
!        once the matrix a has been formed, sqrank should be
!      called once to decompose it.  then for each right hand
!      side, b, sqrlss should be called once to obtain the
!      solution and residual.
!
!  reference(s)
!      dongarra, et al, linpack users guide, siam, 1979
!
  integer lda
  integer n
!
  real a(lda,1)
  real b(*)
  integer info
  integer j
  integer jpvt(1)
  integer k
  integer kr
  integer m
  real qraux(*)
  real rsd(*)
  real t
  real x(n)
!
  if ( kr /= 0 ) then
    call sqrsl(a,lda,m,kr,qraux,b,rsd,rsd,x,rsd,rsd,110,info)
  end if

  jpvt(1:n) = - jpvt(1:n)

  do j = 1, n
    if (j > kr) x(j) = 0.0E+00
  end do

  do j = 1, n

    if (jpvt(j) <= 0) then

      k = -jpvt(j)
      jpvt(j) = k

      do while ( k /= j )
        t = x(j)
        x(j) = x(k)
        x(k) = t
        jpvt(k) = -jpvt(k)
        k = jpvt(k)
      end do

    end if

  end do

  return
end
subroutine sqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
!
!*******************************************************************************
!
!! SQRSL applies the output of sqrdc to compute coordinate
!     transformations, projections, and least squares solutions.
!     for k <= min(n,p), let xk be the matrix
!
!            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
!
!     formed from columnns jpvt(1), ... ,jpvt(k) of the original
!     n x p matrix x that was input to sqrdc (if no pivoting was
!     done, xk consists of the first k columns of x in their
!     original order).  sqrdc produces a factored orthogonal matrix q
!     and an upper triangular matrix r such that
!
!              xk = q * (r)
!                       (0)
!
!     this information is contained in coded form in the arrays
!     x and qraux.
!
!     on entry
!
!        x      real(ldx,p)
!               x contains the output of sqrdc.
!
!        ldx    integer
!               ldx is the leading dimension of the array x.
!
!        n      integer
!               n is the number of rows of the matrix xk.  it must
!               have the same value as n in sqrdc.
!
!        k      integer
!               k is the number of columns of the matrix xk.  k
!               must not be greater than min(n,p), where p is the
!               same as in the calling sequence to sqrdc.
!
!        qraux  real(p)
!               qraux contains the auxiliary output from sqrdc.
!
!        y      real(n)
!               y contains an n-vector that is to be manipulated
!               by sqrsl.
!
!        job    integer
!               job specifies what is to be computed.  job has
!               the decimal expansion abcde, with the following
!               meaning.
!
!                    if a /= 0, compute qy.
!                    if b,c,d, or e /= 0, compute qty.
!                    if c /= 0, compute b.
!                    if d /= 0, compute rsd.
!                    if e /= 0, compute xb.
!
!               note that a request to compute b, rsd, or xb
!               automatically triggers the computation of qty, for
!               which an array must be provided in the calling
!               sequence.
!
!     on return
!
!        qy     real(n).
!               qy contains q*y, if its computation has been
!               requested.
!
!        qty    real(n).
!               qty contains trans(q)*y, if its computation has
!               been requested.  here trans(q) is the
!               transpose of the matrix q.
!
!        b      real(k)
!               b contains the solution of the least squares problem
!
!                    minimize norm2(y - xk*b),
!
!               if its computation has been requested.  (note that
!               if pivoting was requested in sqrdc, the j-th
!               component of b will be associated with column jpvt(j)
!               of the original matrix x that was input into sqrdc.)
!
!        rsd    real(n).
!               rsd contains the least squares residual y - xk*b,
!               if its computation has been requested.  rsd is
!               also the orthogonal projection of y onto the
!               orthogonal complement of the column space of xk.
!
!        xb     real(n).
!               xb contains the least squares approximation xk*b,
!               if its computation has been requested.  xb is also
!               the orthogonal projection of y onto the column space
!               of x.
!
!        info   integer.
!               info is zero unless the computation of b has
!               been requested and r is exactly singular.  in
!               this case, info is the index of the first zero
!               diagonal element of r and b is left unaltered.
!
!     the parameters qy, qty, b, rsd, and xb are not referenced
!     if their computation is not requested and in this case
!     can be replaced by dummy variables in the calling program.
!     to save storage, the user may in some cases use the same
!     array for different parameters in the calling sequence.  a
!     frequently occuring example is when one wishes to compute
!     any of b, rsd, or xb and does not need y or qty.  in this
!     case one may identify y, qty, and one of b, rsd, or xb, while
!     providing separate arrays for anything else that is to be
!     computed.  thus the calling sequence
!
!          call sqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
!
!     will result in the computation of b and rsd, with rsd
!     overwriting y.  more generally, each item in the following
!     list contains groups of permissible identifications for
!     a single callinng sequence.
!
!          1. (y,qty,b) (rsd) (xb) (qy)
!
!          2. (y,qty,rsd) (b) (xb) (qy)
!
!          3. (y,qty,xb) (b) (rsd) (qy)
!
!          4. (y,qy) (qty,b) (rsd) (xb)
!
!          5. (y,qy) (qty,rsd) (b) (xb)
!
!          6. (y,qy) (qty,xb) (b) (rsd)
!
!     in any group the value returned in the array allocated to
!     the group corresponds to the last member of the group.
!
!     linpack.  this version dated 08/14/78 .
!     g. w. stewart, university of maryland, argonne national lab.
!
!  references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
!                 *linpack users  guide*, siam, 1979.
!
  integer k
  integer ldx
!
  real b(k)
  logical cb
  logical cqty
  logical cqy
  logical cr
  logical cxb
  integer i
  integer info
  integer j
  integer jj
  integer job
  integer ju
  integer n
  real qraux(1)
  real qty(1)
  real qy(1)
  real rsd(1)
  real sdot
  real t
  real temp
  real x(ldx,1)
  real xb(1)
  real y(1)
!
!  set info flag.
!
  info = 0
!
!  determine what is to be computed.
!
  cqy = job/10000 /= 0
  cqty = mod(job,10000) /= 0
  cb = mod(job,1000)/100 /= 0
  cr = mod(job,100)/10 /= 0
  cxb = mod(job,10) /= 0
  ju = min(k,n-1)
!
!  special action when n=1.
!
  if (ju == 0) then
     if (cqy) qy(1) = y(1)
     if (cqty) qty(1) = y(1)
     if (cxb) xb(1) = y(1)
     if ( cb ) then
        if (x(1,1) == 0.0) then
           info = 1
        else
           b(1) = y(1)/x(1,1)
        end if
     end if
     if (cr) rsd(1) = 0.0E+00
     go to 250
  end if
!
!  set up to compute qy or qty.
!
  if (cqy) then
    qy(1:n) = y(1:n)
  end if

  if (cqty) then
    qty(1:n) = y(1:n)
  end if
!
!  Compute QY.
!
  if ( cqy ) then

    do jj = 1, ju
      j = ju - jj + 1
      if ( qraux(j) /= 0.0E+00 ) then
        temp = x(j,j)
        x(j,j) = qraux(j)
        t = -sdot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
        call saxpy(n-j+1,t,x(j,j),1,qy(j),1)
        x(j,j) = temp
      end if
    end do

  end if
!
!  Compute Q' * Y.
!
  if ( cqty ) then

        do j = 1, ju
           if ( qraux(j) /= 0.0) then
              temp = x(j,j)
              x(j,j) = qraux(j)
              t = -sdot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
              call saxpy(n-j+1,t,x(j,j),1,qty(j),1)
              x(j,j) = temp
           end if
        end do

  end if
!
!  set up to compute b, rsd, or xb.
!
  if ( cb ) then
    b(1:k) = qty(1:k)
  end if

  if ( cxb ) then
    xb(1:k) = qty(1:k)
  end if

  if ( cr .and. k < n ) then
    rsd(k+1:n) = qty(k+1:n)
  end if

  if ( cxb .and. k+1 <= n ) then
    xb(k+1:n) = 0.0E+00
  end if

  if ( cr ) then
    rsd(1:k) = 0.0E+00
  end if
!
!  compute b.
!
  if ( cb ) then

        do jj = 1, k

           j = k - jj + 1

           if (x(j,j) == 0.0) then
              info = j
              go to 180
           end if

           b(j) = b(j)/x(j,j)
           if (j /= 1) then
              t = -b(j)
              call saxpy(j-1,t,x(1,j),1,b,1)
           end if

        end do

  180       continue

  end if

     if (.not.cr .and. .not.cxb) go to 240
!
!  Compute rsd or xb as required.
!
        do jj = 1, ju

           j = ju - jj + 1
           if (qraux(j) == 0.0E+00 ) go to 220

              temp = x(j,j)
              x(j,j) = qraux(j)

              if ( cr ) then
                 t = -sdot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                 call saxpy(n-j+1,t,x(j,j),1,rsd(j),1)
              end if

              if (.not.cxb) go to 210
                 t = -sdot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                 call saxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
              x(j,j) = temp
  220          continue
        end do

  240    continue
  250 continue

  return
end
subroutine srot ( n, x, incx, y, incy, c, s )
!
!*******************************************************************************
!
!! SROT applies a plane rotation.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input/output, real X(*), one of the vectors to be rotated.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, real Y(*), one of the vectors to be rotated.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
!    Input, real C, S, parameters (presumably the cosine and sine of
!    some angle) that define a plane rotation.
!
  real c
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer n
  real s
  real stemp
  real x(*)
  real y(*)
!
  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    do i = 1, n
      stemp = c * x(i) + s * y(i)
      y(i) = c * y(i) - s * x(i)
      x(i) = stemp
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      stemp = c * x(ix) + s * y(iy)
      y(iy) = c * y(iy) - s * x(ix)
      x(ix) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine srotg ( sa, sb, c, s )
!
!*******************************************************************************
!
!! SROTG constructs a Givens plane rotation.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input/output, real SA, SB, ...
!
!    Output, real C, S, ...
!
  real c
  real r
  real roe
  real s
  real sa
  real sb
  real scale
  real z
!
  if ( abs ( sa ) > abs ( sb ) ) then
    roe = sa
  else
    roe = sb
  end if

  scale = abs ( sa ) + abs ( sb )

  if ( scale == 0.0E+00 ) then
    c = 1.0E+00
    s = 0.0E+00
    r = 0.0E+00
  else
    r = scale * sqrt ( ( sa / scale )**2 + ( sb / scale )**2 )
    r = sign ( 1.0, roe ) * r
    c = sa / r
    s = sb / r
  end if

  if ( abs ( c ) > 0.0E+00 .and. abs ( c ) <= s ) then
    z = 1.0E+00 / c
  else
    z = s
  end if

  sa = r
  sb = z

  return
end
subroutine sscal ( n, sa, x, incx )
!
!*******************************************************************************
!
!! SSCAL scales a vector by a constant.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real SA, the multiplier.
!
!    Input/output, real X(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive entries of X.
!
  integer i
  integer incx
  integer ix
  integer m
  integer n
  real sa
  real x(*)
!
  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    m = mod ( n, 5 )

    x(1:m) = sa * x(1:m)

    do i = m+1, n, 5
      x(i)   = sa * x(i)
      x(i+1) = sa * x(i+1)
      x(i+2) = sa * x(i+2)
      x(i+3) = sa * x(i+3)
      x(i+4) = sa * x(i+4)
    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa * x(ix)
      ix = ix + incx
    end do

  end if

  return
end
subroutine ssvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
!
!*******************************************************************************
!
!! SSVDC computes the singular value decomposition of a rectangular matrix.
!
!
!  The routine reduces a real nxp matrix x by
!     orthogonal transformations u and v to diagonal form.  the
!     diagonal elements s(i) are the singular values of x.  the
!     columns of u are the corresponding left singular vectors,
!     and the columns of v the right singular vectors.
!
!     on entry
!
!         x         real(ldx,p), where ldx >= n.
!                   x contains the matrix whose singular value
!                   decomposition is to be computed.  x is
!                   destroyed by ssvdc.
!
!         ldx       integer
!                   ldx is the leading dimension of the array x.
!
!         n         integer
!                   n is the number of rows of the matrix x.
!
!         p         integer
!                   p is the number of columns of the matrix x.
!
!         ldu       integer
!                   ldu is the leading dimension of the array u.
!                   (see below).
!
!         ldv       integer
!                   ldv is the leading dimension of the array v.
!                   (see below).
!
!         work      real(n)
!                   work is a scratch array.
!
!         job       integer
!                   job controls the computation of the singular
!                   vectors.  it has the decimal expansion ab
!                   with the following meaning
!
!                        a == 0  do not compute the left singular
!                                  vectors.
!                        a == 1  return the n left singular vectors
!                                  in u.
!                        a >= 2  return the first min(n,p) singular
!                                  vectors in u.
!                        b == 0  do not compute the right singular
!                                  vectors.
!                        b == 1  return the right singular vectors
!                                  in v.
!
!     on return
!
!         s         real(mm), where mm=min(n+1,p).
!                   the first min(n,p) entries of s contain the
!                   singular values of x arranged in descending
!                   order of magnitude.
!
!         e         real(p).
!                   e ordinarily contains zeros.  however, see the
!                   discussion of info for exceptions.
!
!         u         real(ldu,k), where ldu >= n.  if joba == 1, then
!                                   k == n.  if joba >= 2 , then
!                                   k == min(n,p).
!                   u contains the matrix of right singular vectors.
!                   u is not referenced if joba == 0.  if n <= p
!                   or if joba == 2, then u may be identified with x
!                   in the subroutine call.
!
!         v         real(ldv,p), where ldv >= p.
!                   v contains the matrix of right singular vectors.
!                   v is not referenced if job == 0.  if p <= n,
!                   then v may be identified with x in the
!                 subroutine call.
!
!         info      integer.
!                   the singular values (and their corresponding
!                   singular vectors) s(info+1),s(info+2),...,s(m)
!                   are correct (here m=min(n,p)).  thus if
!                   info == 0, all the singular values and their
!                   vectors are correct.  in any event, the matrix
!                   b = trans(u)*x*v is the bidiagonal matrix
!                   with the elements of s on its diagonal and the
!                   elements of e on its super-diagonal (trans(u)
!                   is the transpose of u).  thus the singular
!                   values of x and b are the same.
!
!     linpack.  this version dated 03/19/79 .
!     g. w. stewart, university of maryland, argonne national lab.
!
  integer ldu
  integer ldv
  integer ldx
!
  real b
  real c
  real cs
  real e(*)
  real el
  real emm1
  real f
  real g
  integer i
  integer info
  integer iter
  integer j
  integer job
  integer jobu
  integer k
  integer kase
  integer kk
  integer n
  integer p
  real s(*)
  real sdot
  real t
  real t1
  real test
  real u(ldu,*)
  real v(ldv,*)
  logical wantu
  logical wantv
  real work(*)
  real x(ldx,*)
  real ztest

  integer l,ll,lls,lm1,lp1,ls,lu,m,maxit
  integer mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
  real snrm2,scale,shift,sl,sm,sn,smm1
!
!  set the maximum number of iterations.
!
  maxit = 30
!
!  determine what is to be computed.
!
  wantu = .false.
  wantv = .false.
  jobu = mod(job,100)/10
  ncu = n
  if (jobu > 1) ncu = min(n,p)
  if (jobu /= 0) wantu = .true.
  if (mod(job,10) /= 0) wantv = .true.
!
!     reduce x to bidiagonal form, storing the diagonal elements
!     in s and the super-diagonal elements in e.
!
  info = 0
  nct = min(n-1,p)
  nrt = max(0,min(p-2,n))
  lu = max(nct,nrt)

  do l = 1, lu

     lp1 = l + 1

     if (l <= nct) then
!
!           compute the transformation for the l-th column and
!           place the l-th diagonal in s(l).
!
        s(l) = snrm2(n-l+1,x(l,l),1)

        if (s(l) /= 0.0) then
           if (x(l,l) /= 0.0) s(l) = sign(s(l),x(l,l))
           call sscal(n-l+1,1.0/s(l),x(l,l),1)
           x(l,l) = 1.0E+00 + x(l,l)
        end if

        s(l) = -s(l)

     end if

     if (p < lp1) go to 50

     do j = lp1, p
!
!  apply the transformation.
!
        if ( l <= nct ) then
          if ( s(l) /= 0.0E+00 ) then
            t = -sdot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
            call saxpy(n-l+1,t,x(l,l),1,x(l,j),1)
          end if
        end if
!
!           place the l-th row of x into  e for the
!           subsequent calculation of the row transformation.
!
        e(j) = x(l,j)

      end do

   50    continue
!
!  place the transformation in u for subsequent back multiplication.
!
     if ( wantu .and. l <= nct ) then
       u(l:n,l) = x(l:n,l)
     end if

     if (l > nrt) go to 150
!
!  compute the l-th row transformation and place the
!  l-th super-diagonal in e(l).
!
        e(l) = snrm2(p-l,e(lp1),1)

        if (e(l) /= 0.0) then
           if (e(lp1) /= 0.0) e(l) = sign(e(l),e(lp1))
           call sscal(p-l,1.0/e(l),e(lp1),1)
           e(lp1) = 1.0E+00 + e(lp1)
        end if

        e(l) = -e(l)
        if (lp1 > n .or. e(l) == 0.0) go to 120
!
!  apply the transformation.
!
           work(lp1:n) = 0.0E+00

           do j = lp1, p
             call saxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
           end do

           do j = lp1, p
             call saxpy(n-l,-e(j)/e(lp1),work(lp1),1,x(lp1,j),1)
           end do

  120       continue
!
!  Place the transformation in v for subsequent back multiplication.
!
         if ( wantv ) then
           v(lp1:p,l) = e(lp1:p)
         end if

  150    continue

  end do
!
!  Set up the final bidiagonal matrix or order m.
!
  m = min(p,n+1)
  nctp1 = nct + 1
  nrtp1 = nrt + 1
  if (nct < p) s(nctp1) = x(nctp1,nctp1)
  if (n < m) s(m) = 0.0E+00
  if (nrtp1 < m) e(nrtp1) = x(nrtp1,m)
  e(m) = 0.0E+00
!
!  Generate u.
!
  if ( wantu ) then

     do j = nctp1, ncu
       u(1:n,j) = 0.0E+00
       u(j,j) = 1.0E+00
     end do

     do ll = 1, nct

       l = nct - ll + 1

       if (s(l) /= 0.0) then

           lp1 = l + 1

           do j = l+1, ncu
             t = -sdot(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
             call saxpy(n-l+1,t,u(l,l),1,u(l,j),1)
           end do

           call sscal ( n-l+1, -1.0, u(l,l), 1 )
           u(l,l) = 1.0E+00 + u(l,l)
           u(1:l-1,l) = 0.0E+00

        else

           u(1:n,l) = 0.0E+00
           u(l,l) = 1.0E+00

        end if

      end do

  end if
!
!  Generate v.
!
  if ( wantv ) then

    do ll = 1, p

      l = p - ll + 1
      lp1 = l + 1

      if ( l <= nrt ) then

        if ( e(l) /= 0.0E+00 ) then

          do j = lp1, p
            t = -sdot(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
            call saxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
          end do

        end if

      end if

      v(1:p,l) = 0.0E+00
      v(l,l) = 1.0E+00

    end do

  end if
!
!  main iteration loop for the singular values.
!
  mm = m
  iter = 0

  360 continue
!
!  quit if all the singular values have been found.
!
     if (m == 0) go to 620
!
!  if too many iterations have been performed, set flag and return.
!
     if ( iter >= maxit ) then
       info = m
       go to 620
     end if
!
!  this section of the program inspects for
!  negligible elements in the s and e arrays.  on
!  completion the variables kase and l are set as follows.
!
!           kase = 1     if s(m) and e(l-1) are negligible and l<m
!           kase = 2     if s(l) is negligible and l<m
!           kase = 3     if e(l-1) is negligible, l<m, and
!                        s(l), ..., s(m) are not negligible (qr step).
!           kase = 4     if e(m-1) is negligible (convergence).
!
     do ll = 1, m

        l = m - ll
        if (l == 0) go to 400
        test = abs(s(l)) + abs(s(l+1))
        ztest = test + abs(e(l))

        if (ztest == test) then
           e(l) = 0.0E+00
           go to 400
        end if

     end do

  400    continue
     if (l /= m - 1) go to 410
        kase = 4
     go to 480
  410    continue
        lp1 = l + 1
        mp1 = m + 1

        do lls = lp1, mp1
           ls = m - lls + lp1
           if (ls == l) go to 440
           test = 0.0E+00
           if (ls /= m) test = test + abs(e(ls))
           if (ls /= l + 1) test = test + abs(e(ls-1))
           ztest = test + abs(s(ls))
           if (ztest /= test) go to 420
              s(ls) = 0.0E+00
              go to 440
  420          continue
        end do

  440       continue
        if (ls /= l) go to 450
           kase = 3
        go to 470
  450       continue
        if (ls /= m) go to 460
           kase = 1
        go to 470
  460       continue
           kase = 2
           l = ls
  470       continue
  480    continue
     l = l + 1
!
!        perform the task indicated by kase.
!
     go to (490,520,540,570), kase
!
!  deflate negligible s(m).
!
  490    continue
        mm1 = m - 1
        f = e(m-1)
        e(m-1) = 0.0E+00
        do kk = l, mm1
           k = mm1 - kk + l
           t1 = s(k)
           call srotg(t1,f,cs,sn)
           s(k) = t1
           if (k == l) go to 500
              f = -sn*e(k-1)
              e(k-1) = cs*e(k-1)
  500          continue
           if (wantv) call srot(p,v(1,k),1,v(1,m),1,cs,sn)
        end do

     go to 610
!
!  Split at negligible s(l).
!
  520    continue
        f = e(l-1)
        e(l-1) = 0.0E+00

        do k = l, m
           t1 = s(k)
           call srotg(t1,f,cs,sn)
           s(k) = t1
           f = -sn*e(k)
           e(k) = cs*e(k)
           if (wantu) call srot(n,u(1,k),1,u(1,l-1),1,cs,sn)
        end do

     go to 610
!
!  perform one qr step.
!
  540    continue
!
!  calculate the shift.
!
        scale = max ( abs(s(m)),abs(s(m-1)),abs(e(m-1)),abs(s(l)),abs(e(l)) )

        sm = s(m)/scale
        smm1 = s(m-1)/scale
        emm1 = e(m-1)/scale
        sl = s(l)/scale
        el = e(l)/scale
        b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0E+00
        c = (sm*emm1)**2
        shift = 0.0E+00

        if ( b /= 0.0E+00 .or. c /= 0.0E+00 ) then
           shift = sqrt(b**2+c)
           if (b < 0.0) shift = -shift
           shift = c/(b + shift)
        end if

        f = (sl + sm)*(sl - sm) - shift
        g = sl*el
!
!  chase zeros.
!
        mm1 = m - 1

        do k = l, m-1
           call srotg(f,g,cs,sn)
           if (k /= l) e(k-1) = f
           f = cs*s(k) + sn*e(k)
           e(k) = cs*e(k) - sn*s(k)
           g = sn*s(k+1)
           s(k+1) = cs*s(k+1)
           if (wantv) call srot(p,v(1,k),1,v(1,k+1),1,cs,sn)
           call srotg(f,g,cs,sn)
           s(k) = f
           f = cs*e(k) + sn*s(k+1)
           s(k+1) = -sn*e(k) + cs*s(k+1)
           g = sn*e(k+1)
           e(k+1) = cs*e(k+1)
           if (wantu .and. k < n) then
             call srot(n,u(1,k),1,u(1,k+1),1,cs,sn)
           end if
        end do

        e(m-1) = f
        iter = iter + 1
     go to 610
!
!  convergence.
!
  570    continue
!
!  make the singular value  positive.
!
        if (s(l) >= 0.0) go to 580
           s(l) = -s(l)
           if (wantv) call sscal(p,-1.0,v(1,l),1)
  580       continue
!
!  order the singular value.
!
  590       if (l == mm) go to 600
           if (s(l) >= s(l+1)) go to 600
           t = s(l)
           s(l) = s(l+1)
           s(l+1) = t
           if (wantv .and. l < p) then
             call sswap(p,v(1,l),1,v(1,l+1),1)
           end if

           if (wantu .and. l < n) then
             call sswap(n,u(1,l),1,u(1,l+1),1)
           end if

           l = l + 1
        go to 590
  600       continue
        iter = 0
        m = m - 1
  610    continue
  go to 360
  620 continue
  return
end
subroutine sswap ( n, x, incx, y, incy )
!
!*******************************************************************************
!
!! SSWAP interchanges two vectors.
!
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input/output, real X(*), one of the vectors to swap.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Input/output, real Y(*), one of the vectors to swap.
!
!    Input, integer INCY, the increment between successive elements of Y.
!
  integer i
  integer incx
  integer incy
  integer ix
  integer iy
  integer m
  integer n
  real stemp
  real x(*)
  real y(*)
!
  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      stemp = x(i)
      x(i) = y(i)
      y(i) = stemp
    end do

    do i = m+1, n, 3

      stemp = x(i)
      x(i) = y(i)
      y(i) = stemp

      stemp = x(i + 1)
      x(i + 1) = y(i + 1)
      y(i + 1) = stemp

      stemp = x(i + 2)
      x(i + 2) = y(i + 2)
      y(i + 2) = stemp

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      stemp = x(ix)
      x(ix) = y(iy)
      y(iy) = stemp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine tregup(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl, &
  dlt,iretcd,xplsp,fplsp,xpls,fpls,mxtake,ipr,method,udiag)
!
!*******************************************************************************
!
!! TREGUP decides whether to accept the next optimization iterate.
!
!
!  Parameters:
!
! nr           --> row dimension of matrix
! n            --> dimension of problem
! x(n)         --> old iterate x[k-1]
! f            --> function value at old iterate, f(x)
! g(n)         --> gradient at old iterate, g(x), or approximate
! a(n,n)       --> cholesky decomposition of hessian in
!                  lower triangular part and diagonal.
!                  hessian or approx in upper triangular part
! fcn          --> name of subroutine to evaluate function
! sc(n)        --> current step
! sx(n)        --> diagonal scaling matrix for x
! nwtake       --> boolean, =.true. if newton step taken
! stepmx       --> maximum allowable step size
! steptl       --> relative step size at which successive iterates
!                  considered close enough to terminate algorithm
! dlt         <--> trust region radius
! iretcd      <--> return code
!                    =0 xpls accepted as next iterate;
!                       dlt trust region for next iteration.
!                    =1 xpls unsatisfactory but accepted as next iterate
!                       because xpls-x < smallest allowable
!                       step length.
!                    =2 f(xpls) too large.  continue current iteration
!                       with new reduced dlt.
!                    =3 f(xpls) sufficiently small, but quadratic model
!                       predicts f(xpls) sufficiently well to continue
!                       current iteration with new doubled dlt.
! xplsp(n)    <--> workspace [value needs to be retained between
!                  succesive calls of k-th global step]
! fplsp       <--> [retain value between successive calls]
! xpls(n)     <--  new iterate x[k]
! fpls        <--function value at new iterate, f(xpls)
! mxtake      <--  boolean flag indicating step of maximum length used
! ipr          --> device to which to send output
! method       --> algorithm to use to solve minimization problem
!                    =1 line search
!                    =2 double dogleg
!                    =3 more-hebdon
! udiag(n)     --> diagonal of hessian in a(.,.)
!
  integer n
  integer nr
!
  real a(nr,n)
  real dlt
  real dltf
  real dltfp
  real dltmp
  real dltp
  real f
  real fpls
  real fplsp
  real g(n)
  integer i
  integer ipr
  integer iretcd
  integer j
  integer method
  logical mxtake
  logical nwtake
  real rln
  real sc(n)
  real sdot
  real slp
  real stepmx
  real steptl
  real sx(n)
  real temp
  real udiag(n)
  real x(n)
  real xpls(n)
  real xplsp(n)
!
  external fcn
!
  mxtake = .false.
  xpls(1:n) = x(1:n) + sc(1:n)
  call fcn ( n, xpls, fpls )
  dltf = fpls - f
  slp = sdot ( n, g, 1, sc, 1 )
!
! next statement added for case of compilers which do not optimize
! evaluation of next "if" statement (in which case fplsp could be
! undefined).
!
  if ( iretcd==4) fplsp=0.0E+00

  if ( iretcd/=3 .or. (fpls<fplsp .and. dltf<= 1.e-4*slp)) then
    go to 130
  end if
!     if ( iretcd==3 .and. (fpls>=fplsp .or. dltf> 1.e-4*slp))
!     then
!
!       reset xpls to xplsp and terminate global step
!
    iretcd=0
    xpls(1:n)=xplsp(1:n)
    fpls=fplsp
    dlt=.5*dlt
    go to 230
!     else
!
!  fpls too large
!
  130   if ( dltf<= 1.e-4*slp) go to 170
!       if ( dltf> 1.e-4*slp)
!       then
      rln=0.
      do i=1,n
        rln=max(rln,abs(sc(i))/max(abs(xpls(i)),1./sx(i)))
      end do

      if ( rln>=steptl) go to 150
!         if ( rln<steptl)
!         then
!
!  cannot find satisfactory xpls sufficiently distinct from x
!
        iretcd=1
        go to 230
!         else
!
!  reduce trust region and continue global step
!
  150       iretcd=2
        dltmp=-slp*dlt/(2.*(dltf-slp))
        if ( dltmp>= .1*dlt) go to 155
!           if ( dltmp< .1*dlt)
!           then
          dlt=.1*dlt
          go to 160
!           else
  155         dlt=dltmp
!           end if
  160       continue
        go to 230
!         end if
!       else
!
!  fpls sufficiently small
!
  170     continue
      dltfp=0.
      if (method == 3) go to 180
!
!         if (method == 2)
!         then
!
      do i = 1, n
         temp = 0.0E+00
         do j = i, n
            temp = temp + (a(j, i)*sc(j))
         end do
         dltfp = dltfp + temp*temp
      end do

      go to 190

  180     continue

      do i = 1, n
         dltfp = dltfp + udiag(i)*sc(i)*sc(i)

         temp = 0.0E+00
         do j = i+1, n
            temp = temp + a(i, j)*sc(i)*sc(j)
         end do
         dltfp = dltfp + 2.0*temp
      end do

  190     dltfp = slp + dltfp/2.0E+00
      if ( iretcd==2 .or. (abs(dltfp-dltf)> .1*abs(dltf)) &
         .or. nwtake .or. (dlt> .99*stepmx)) go to 210
!         if ( iretcd/=2 .and. (abs(dltfp-dltf) <= .1*abs(dltf))
!    +         .and. (.not.nwtake) .and. (dlt<= .99*stepmx))
!         then
!
!           double trust region and continue global step
!
        iretcd=3
        xplsp(1:n)=xpls(1:n)
        fplsp=fpls
        dlt=min(2.*dlt,stepmx)
        go to 230
!         else
!
!           accept xpls as next iterate.  choose new trust region.
!
  210       continue
        iretcd=0
        if ( dlt> .99*stepmx) mxtake=.true.
        if ( dltf< .1*dltfp) go to 220
!           if ( dltf>= .1*dltfp)
!           then
!
!  decrease trust region for next iteration
!
          dlt=.5*dlt
          go to 230
!           else
!
!             check whether to increase trust region for next iteration
!
  220         if ( dltf<= .75*dltfp) dlt=min(2.*dlt,stepmx)
!           end if
!         end if
!       end if
!     end if
  230 continue
  return
end
subroutine uncmin (n,x0,fcn,x,f,info,w,lw)
!
!*******************************************************************************
!
!! UNCMIN minimizes a smooth nonlinear function of n variables.
!
!
!     a subroutine that computes the function value at any point
!     must be supplied.  derivative values are not required.
!     this subroutine provides the simplest interface to the uncmin
!     minimization package.  user has no control over options.
!
!  Discussion:
!
!     this routine uses a quasi-newton algorithm with line search
!     to minimize the function represented by the subroutine fcn.
!     at each iteration, the nonlinear function is approximated
!     by a quadratic function derived from a taylor series.
!     the quadratic function is minimized to obtain a search direction,
!     and an approximate minimum of the nonlinear function along
!     the search direction is found using a line search.  the
!     algorithm computes an approximation to the second derivative
!     matrix of the nonlinear function using quasi-newton techniques.
!
!     the uncmin package is quite general, and provides many options
!     for the user.  however, this subroutine is designed to be
!     easy to use, with few choices allowed.  for example:
!
!     1.  only function values need be computed.  first derivative
!     values are obtained by finite-differencing.  this can be
!     very costly when the number of variables is large.
!
!     2.  it is assumed that the function values can be obtained
!     accurately (to an accuracy comparable to the precision of
!     the computer arithmetic).
!
!     3.  at most 150 iterations are allowed.
!
!     4.  it is assumed that the function values are well-scaled,
!     that is, that the optimal function value is not pathologically
!     large or small.
!
!     for more information, see the reference listed below.
!
!  Parameters:
!
! n            --> integer
!                  dimension of problem
! x0(n)        --> real
!                  initial estimate of minimum
! fcn          --> name of routine to evaluate minimization function.
!                  must be declared external in calling routine, and
!                  have calling sequence
!                    subroutine fcn(n, x, f)
!                  with n and x as here, f the computed function value.
! x(n)        <--  real
!                  local minimum
! f           <--  real
!                function value at local minimum x
! info        <--  integer
!                  termination code
!                      info =  0:  optimal solution found
!                      info =  1:  terminated with gradient small,
!                                  x is probably optimal
!                      info =  2:  terminated with stepsize small,
!                                  x is probably optimal
!                      info =  3:  lower point cannot be found,
!                                  x is probably optimal
!                      info =  4:  iteration limit (150) exceeded
!                      info =  5:  too many large steps,
!                                function may be unbounded
!                      info = -1:  insufficient workspace
! w(lw)        --> real
!                  workspace
! lw           --> integer
!                  size of workspace, at least n*(n+10)
!
!  references  r.b. schnabel, j.e. koontz, and be.e. weiss, a modular
!                 system of algorithms for unconstrained minimization,
!                 report cu-cs-240-82, comp. sci. dept., univ. of
!                 colorado at boulder, 1982.
!
  integer lw
  integer n
!
  real dlt
  real epsm
  character errmsg*68
  real f
  real fscale
  real gradtl
  integer ia
  integer iagflg
  integer iahflg
  integer iexp
  integer ig
  integer info
  integer ipr
  integer it
  integer itnlim
  integer iw1
  integer iw2
  integer iw3
  integer iw4
  integer iw5
  integer iw6
  integer iw7
  integer iw8
  integer lwmin
  integer method
  integer msg
  integer ndigit
  integer nr
  real stepmx
  real steptl
  real w(lw)
  real x(n)
  real x0(n)
!
  external  fcn, d1fcn, d2fcn
!
! subdivide workspace
!
  ig  = 1
  it  = ig  + n
  iw1 = it  + n
  iw2 = iw1 + n
  iw3 = iw2 + n
  iw4 = iw3 + n
  iw5 = iw4 + n
  iw6 = iw5 + n
  iw7 = iw6 + n
  iw8 = iw7 + n
  ia  = iw8 + n
  lwmin = ia + n*n-1

  if (lwmin > lw) then
    info = -1
    write ( *, * ) ' '
    write ( *, * ) 'UNCMIN - Fatal error!'
    write ( *, * ) '  Insufficient workspace.'
    write ( *, * ) '  LW < LWMIN.'
    write ( *, * ) '  LW = ', lw
    write ( *, * ) '  LWMIN = ', lwmin
    stop
  end if
!
!  set up parameters for optdrv
!
!  parameters that should not be changed when using condensed code
!
! nr     = parameter used to divide workspace
! method = 1 (line search) -- do not change
! msg    = 9 => no printing, n=1 allowed
! iahflg = 1 => analytic hessian  supplied (0 otherwise)
! ipr    = device for output (irrelevant in current version)
! dlt    = (irrelevant parameter for method = 1)
! epsm   = machine epsilon
!
  nr = n
  method = 1
  msg = 9
  iahflg = 0
  ipr = 6
  dlt = -1.0E+00
  epsm = epsilon ( epsm )
!
! parameters that may be changed:
!
! iexp   = 1 => function expensive to evaluate (iexp = 0 => cheap)
! iagflg = 1 => analytic gradient supplied (0 otherwise)
! ndigit = -1 => optdrv assumes f is fully accurate
! itnlim = 150 = maximum number of iterations allowed
! gradtl = zero tolerance for gradient, for convergence tests
! stepmx = maximum allowable step size
! steptl = zero tolerance for step, for convergence tests
! fscale = typical order-of-magnitude size of function
! typsiz = typical order-of-magnitude size of x (stored in w(lt))
!
  iexp = 1
  iagflg = 0
  ndigit = -1
  itnlim = 150
  gradtl = epsm**(1.0/3.0)
  stepmx = 0.0E+00
  steptl = sqrt(epsm)
  fscale = 1.0E+00
  w(it:it+n-1) = 1.0E+00
!
! minimize function
!
  call optdrv (nr, n, x0, fcn, d1fcn, d2fcn, w(it), fscale, &
                  method, iexp, msg, ndigit, itnlim, iagflg, iahflg, &
                  ipr, dlt, gradtl, stepmx, steptl, &
                  x, f, w(ig), info, w(ia), &
                  w(iw1), w(iw2), w(iw3), w(iw4), &
                  w(iw5), w(iw6), w(iw7), w(iw8))

  if (info == 1) then
    write ( *, * ) ' '
    write ( *, * ) 'UNCMIN - Note!'
    write ( *, * ) '  INFO = 1.'
    write ( *, * ) '  The iteration probably converged.'
    write ( *, * ) '  The gradient is very small.'
    return
  end if

  if (info == 2) then
    write ( *, * ) ' '
    write ( *, * ) 'UNCMIN - Note!'
    write ( *, * ) '  INFO = 2.'
    write ( *, * ) '  The iteration probably converged.'
    write ( *, * ) '  The stepsize is very small.'
    return
  end if

  if (info == 3) then
    write ( *, * ) ' '
    write ( *, * ) 'UNCMIN - Warning!'
    write ( *, * ) '  INFO = 3.'
    write ( *, * ) '  Cannot find a point with lower value.'
    write ( *, * ) '  (But not completely happy with the current value.)'
    return
  end if

  if (info == 4) then
    write ( *, * ) ' '
    write ( *, * ) 'UNCMIN - Warning!'
    write ( *, * ) '  INFO = 4.'
    write ( *, * ) '  Too many iterations.'
    return
  end if

  if (info == 5) then
    write ( *, * ) ' '
    write ( *, * ) 'UNCMIN - Warning!'
    write ( *, * ) '  INFO = 5.'
    write ( *, * ) '  Too many large steps.'
    write ( *, * ) '  The function may be unbounded.'
    return
  end if

  return
end
function uni ( )
!
!*******************************************************************************
!
!! UNI generates real uniform random numbers on [0,1).
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!       usage:
!              to initialize the generator
!                   useed = ustart(iseed)
!               where: iseed is any nonzero integer
!                  will return floating point value of iseed.
!
!               subsequently
!                       u = uni()
!                  will return a real uniform on [0,1)
!
!                one initialization is necessary, but any number of evaluations
!                  of  uni in any order, are allowed.
!
!           note: depending upon the value of k (see below), the output
!                       of uni may differ from one machine to another.
!
!           typical usage:
!
!               real u,uni,ustart,useed
!               integer iseed
!c                 set seed
!               iseed = 305
!               useed = ustart(iseed)
!               do i = 1,1000
!                   u = uni()
!               end do
!c                 note: if k=24 (the default, see below) the output value of
!c                           u will be 0.1570390462475...
!               write(*,*) u
!             end
!
!          note on portability: users can choose to run uni in its default
!               mode (requiring no user action) which will generate the same
!               sequence of numbers on any computer supporting floating point
!               numbers with at least 24 bit mantissas, or in a mode that
!               will generate numbers with a longer period on computers with
!               larger mantissas.
!          to exercise this option:  b e f o r e  invoking ustart insert
!               the instruction        ubits = unib(k)      k >= 24
!               where k is the number of bits in the mantissa of your floating
!               point word (k=48 for cray, cyber 205). unib returns the
!               floating point value of k that it actually used.
!                    k input as <= 24, then ubits=24.
!                    k input as > 24, then ubits=float(k)
!               if k>24 the sequence of numbers generated by uni may differ
!               from one computer to another.
!
!
!
!  references  marsaglia g., "comments on the perfect uniform random
!                 number generator", unpublished notes, wash s. u.
!
  real, save :: c = 362436.0E+00 / 16777216.0E+00
  real, parameter :: cd = 7654321.0E+00 / 16777216.0E+00
  real, parameter :: cm = 16777213.0E+00 / 16777216.0E+00
  real, parameter :: csave = 362436.0E+00 / 16777216.0E+00
  integer, save :: i = 17
  integer i1
  integer ii
  integer iseed
  integer, save :: j = 5
  integer j1
  integer jj
  integer, save :: k = 24
  integer k1
  integer kk
  integer l1
  integer m1
  real s
  real t
  real, save, dimension ( 17 ) :: u = (/ &
     0.8668672834288,  0.3697986366357,  0.8008968294805, &
     0.4173889774680,  0.8254561579836,  0.9640965269077, &
     0.4508667414265,  0.6451309529668,  0.1645456024730, &
     0.2787901807898,  0.06761531340295, 0.9663226330820, &
     0.01963343943798, 0.02947398211399, 0.1636231515294, &
     0.3976343250467,  0.2631008574685 /)
  real uni
  real unib
  real ustart
!
!      load data array in case user forgets to initialize.
!      this array is the result of calling uni 100000 times
!         with iseed=305 and k=64.
!
!   basic generator is fibonacci
!
  uni = u(i)-u(j)
  if ( uni<0.0)uni = uni+1.0E+00
  u(i) = uni
  i = i-1
  if ( i==0)i = 17
  j = j-1
  if ( j==0)j = 17
!
!   second generator is congruential
!
  c = c-cd
  if ( c<0.0) c=c+cm
!
!   combination generator
!
  uni = uni-c
  if ( uni<0.0)uni = uni+1.0E+00
  return
!
entry ustart ( iseed )
!
!          set up ...
!          convert iseed to four smallish positive integers.
!
    i1 = mod(abs(iseed),177)+1
    j1 = mod(abs(iseed),167)+1
    k1 = mod(abs(iseed),157)+1
    l1 = mod(abs(iseed),147)+1
!
!              generate random bit pattern in array based on given seed.
!
    do ii = 1,17

      s = 0.0E+00
      t = 0.5
!             do for each of the bits of mantissa of word
!             loop  over k bits, where k is defaulted to 24 but can
!               be changed by user call to unib(k)
!
      do jj = 1,k
        m1 = mod(mod(i1*j1,179)*k1,179)
        i1 = j1
        j1 = k1
        k1 = m1
        l1 = mod(53*l1+1,169)
        if ( mod(l1*m1,64)>=32)s=s+t
        t = .5*t
       end do

       u(ii) = s

    end do

    ustart = real(iseed)
    return
!
entry unib ( kk )

  if ( kk <= 24 ) then
    k = 24
  else
    k = kk
  end if

  unib = real ( k )

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
subroutine xerabt(messg,nmessg)
!
!*******************************************************************************
!
!! XERABT aborts program execution and prints error message.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        ***note*** machine dependent routine
!        xerabt aborts the execution of the program.
!        the error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     description of parameters
!        messg and nmessg are as in xerror, except that nmessg may
!        be zero, in which case no message is being supplied.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  19 mar 1980
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  character*(*) messg
  integer nmessg
!
  stop
end
subroutine xerclr
!
!*******************************************************************************
!
!! XERCLR resets current error number to zero.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        this routine simply resets the current error number to zero.
!        this may be necessary to do in order to determine that
!        a certain error has occurred again since the last time
!        numxer was referenced.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer j4save
  integer junk
!
  junk = j4save ( 1, 0, .true. )

  return
end
subroutine xerctl(messg1,nmessg,nerr,level,kontrl)
!
!*******************************************************************************
!
!! XERCTL allows user control over handling of individual errors.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        allows user control over handling of individual errors.
!        just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to xerctl.
!        if the user has provided his own version of xerctl, he
!        can then override the value of kontrol used in processing
!        this message by redefining its value.
!        kontrl may be set to any value from -2 to 2.
!        the meanings for kontrl are the same as in xsetf, except
!        that the value of kontrl changes only for this message.
!        if kontrl is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     description of parameters
!
!      --input--
!        messg1 - the first word (only) of the error message.
!        nmessg - same as in the call to xerror or xerrwv.
!        nerr   - same as in the call to xerror or xerrwv.
!        level  - same as in the call to xerror or xerrwv.
!        kontrl - the current value of the control flag as set
!                 by a call to xsetf.
!
!      --output--
!        kontrl - the new value of kontrl.  if kontrl is not
!                 defined, it will remain at its original value.
!                 this changed value of control affects only
!                 the current occurrence of the current message.
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer kontrl
  integer level
  character*20 messg1
  integer nerr
  integer nmessg
!
  return
end
subroutine xerdmp
!
!*******************************************************************************
!
!! XERDMP prints the error tables and then clears them.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xerdmp prints the error tables, then clears them.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer kount
!
  call xersav(' ',0,0,0,kount)

  return
end
subroutine xermax ( maxnum )
!
!*******************************************************************************
!
!! XERMAX sets maximum number of times any error message is to be printed.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xermax sets the maximum number of times any message
!        is to be printed.  that is, non-fatal messages are
!        not to be printed after they have occured maxnum times.
!        such non-fatal messages may be printed less than
!        maxnum times even if they occur maxnum times, if error
!        suppression mode (kontrl=0) is ever in effect.
!
!     description of parameter
!      --input--
!        maxnum - the maximum number of times any one message
!              is to be printed.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer j4save
  integer junk
  integer maxnum
!
  junk = j4save ( 4, maxnum, .true. )

  return
end
subroutine xerprt(messg,nmessg)
!
!*******************************************************************************
!
!! XERPRT prints a message on each file indicated by xgetua.
!
!
!     latest revision ---  16 sept 1987
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer ichar
  integer iunit
  integer kunit
  integer last
  integer lenmes
  integer lun(5)
  character*(*) messg
  integer nmessg
  integer nunit
!
!     obtain unit numbers and write line to each unit
!
  call xgetua(lun,nunit)
  lenmes = len(messg)

  do kunit=1,nunit

     iunit = lun(kunit)

     do ichar=1,lenmes,72
        last = min(ichar+71 , lenmes)
        if ( iunit==0 ) then
          write (*,'(1x,a)') messg(ichar:last)
        else
          write (iunit,'(1x,a)') messg(ichar:last)
        end if
    end do

  end do

  return
end
subroutine xerror(messg,nmessg,nerr,level)
!
!*******************************************************************************
!
!! XERROR processes an error (diagnostic) message.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xerror processes a diagnostic message, in a manner
!        determined by the value of level and the current value
!        of the library error control flag, kontrl.
!        (see subroutine xsetf for details.)
!
!     description of parameters
!      --input--
!        messg - the hollerith message to be processed, containing
!                no more than 72 characters.
!        nmessg- the actual number of characters in messg.
!        nerr  - the error number associated with this message.
!                nerr must not be zero.
!        level - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (i.e., it is
!                   non-fatal if xsetf has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!
!     examples
!        call xerror('smooth -- num was zero.',23,1,2)
!        call xerror('integ  -- less than full accuracy achieved.',
!                    43,2,1)
!        call xerror('rooter -- actual zero of f found before interval f
!    1ully collapsed.',65,3,0)
!        call xerror('exp    -- underflows being set to zero.',39,1,-1)
!
!     latest revision ---  19 mar 1980
!     written by ron jones, with slatec common math library subcommittee
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer level
  character*(*) messg
  integer nerr
  integer nmessg
!
  call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)

  return
end
subroutine xerrwv(messg,nmessg,nerr,level,ni,i1,i2,nr,r1,r2)
!
!*******************************************************************************
!
!! XERRWV processes error message allowing 2 integer and two real
!            values to be included in the message.
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xerrwv processes a diagnostic message, in a manner
!        determined by the value of level and the current value
!        of the library error control flag, kontrl.
!        (see subroutine xsetf for details.)
!        in addition, up to two integer values and two real
!        values may be printed along with the message.
!
!     description of parameters
!      --input--
!        messg - the hollerith message to be processed.
!        nmessg- the actual number of characters in messg.
!        nerr  - the error number associated with this message.
!                nerr must not be zero.
!        level - error category.
!                =2 means this is an unconditionally fatal error.
!                =1 means this is a recoverable error.  (i.e., it is
!                   non-fatal if xsetf has been appropriately called.)
!                =0 means this is a warning message only.
!                =-1 means this is a warning message which is to be
!                   printed at most once, regardless of how many
!                   times this call is executed.
!        ni    - number of integer values to be printed. (0 to 2)
!        i1    - first integer value.
!        i2    - second integer value.
!        nr    - number of real values to be printed. (0 to 2)
!        r1    - first real value.
!        r2    - second real value.
!
!     examples
!        call xerrwv('smooth -- num (=i1) was zero.',29,1,2,
!    1   1,num,0,0,0.,0.)
!        call xerrwv('quadxy -- requested error (r1) less than minimum (
!    1r2).,54,77,1,0,0,0,2,errreq,errmin)
!
!     latest revision ---  16 sept 1987
!     written by ron jones, with slatec common math library subcommittee
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  character*37 form
  integer i
  integer i1
  integer i1mach
  integer i2
  integer ifatal
  integer isizei
  integer isizef
  integer iunit
  integer j4save
  integer junk
  integer kdummy
  integer kount
  integer kunit
  integer lerr
  integer level
  character*20 lfirst
  integer lkntrl
  integer llevel
  integer lmessg
  integer lun(5)
  integer maxmes
  character*(*) messg
  integer mkntrl
  integer nerr
  integer ni
  integer nmessg
  integer nr
  integer nunit
  real r1
  real r2
!
!     get flags
!
  lkntrl = j4save(2,0,.false.)
  maxmes = j4save(4,0,.false.)
!     check for valid input
  if ((nmessg>0).and.(nerr/=0).and.(level>=(-1)).and.(level<=2)) then
    go to 10
  end if

     if (lkntrl>0) call xerprt('fatal error in...',17)
     call xerprt('xerror -- invalid input',23)

     if (lkntrl>0) then
       call xerprt('job abort due to fatal error.',29)
     end if
     if (lkntrl>0) call xersav(' ',0,0,0,kdummy)
     call xerabt('xerror -- invalid input',23)
     return
   10 continue
!     record message
  junk = j4save(1,nerr,.true.)
  call xersav(messg,nmessg,nerr,level,kount)
!     let user override
  lfirst = messg
  lmessg = nmessg
  lerr = nerr
  llevel = level
  call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
!
!  Reset to original values.
!
  lmessg = nmessg
  lerr = nerr
  llevel = level
  lkntrl = max(-2,min(2,lkntrl))
  mkntrl = abs(lkntrl)
!     decide whether to print message
  if ((llevel<2).and.(lkntrl==0)) go to 100

  if (((llevel==(-1)).and.(kount>min(1,maxmes))) &
    .or.((llevel==0)   .and.(kount>maxmes)) &
    .or.((llevel==1)   .and.(kount>maxmes).and.(mkntrl==1)) &
    .or.((llevel==2)   .and.(kount>max(1,maxmes)))) then
    go to 100
  end if

     if (lkntrl<=0) go to 20

        call xerprt(' ',1)

        if (llevel==(-1)) then
          call xerprt &
      ('warning message...this message will only be printed once.',57)
        end if

        if (llevel==0) call xerprt('warning in...',13)

        if (llevel==1) call xerprt('recoverable error in...',23)

        if (llevel==2) call xerprt('fatal error in...',17)
   20    continue
!        message
     call xerprt(messg,lmessg)
     call xgetua(lun,nunit)
     isizei = log10(real(i1mach(9))) + 1.0E+00
     isizef = log10(real(i1mach(10))**i1mach(11)) + 1.0E+00

     do kunit=1,nunit

        iunit = lun(kunit)

        do i=1,min(ni,2)
           write (form,21) i,isizei
   21          format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
           if ( iunit==0 ) then
             if (i==1) write (*,form) i1
             if (i==2) write (*,form) i2
           else
             if (i==1) write (iunit,form) i1
             if (i==2) write (iunit,form) i2
           end if
        end do

        do i=1,min(nr,2)
           write (form,23) i,isizef+10,isizef
   23          format ('(11x,21hin above message, r',i1,'=,e',i2,'.',i2,')')
           if ( iunit==0 ) then
             if (i==1) write (*,form) r1
             if (i==2) write (*,form) r2
           else
             if (i==1) write (iunit,form) r1
             if (i==2) write (iunit,form) r2
           end if
        end do

        if (lkntrl<=0) go to 40
!              error number
           if ( iunit==0 ) then
             write(*,30) lerr
           else
             write (iunit,30) lerr
           end if
   30          format (15h error number =,i10)
   40       continue

     end do
!
!        trace-back

  100 continue
  ifatal = 0
  if ((llevel==2).or.((llevel==1).and.(mkntrl==2))) then
    ifatal = 1
  end if

!     quit here if message is not fatal
  if (ifatal<=0) return
  if ((lkntrl<=0).or.(kount>max(1,maxmes))) go to 120
!        print reason for abort
     if (llevel==1) call xerprt &
        ('job abort due to unrecovered error.',35)

     if (llevel==2) then
       call xerprt('job abort due to fatal error.',29)
     end if

!        print error summary
     call xersav(' ',-1,0,0,kdummy)
  120 continue
!     abort
  if ((llevel==2).and.(kount>max(1,maxmes))) lmessg = 0
  call xerabt(messg,lmessg)
  return
end
subroutine xersav(messg,nmessg,nerr,level,icount)
!
!*******************************************************************************
!
!! XERSAV records that an error occurred.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        record that this error occurred.
!
!     description of parameters
!     --input--
!       messg, nmessg, nerr, level are as in xerror,
!       except that when nmessg=0 the tables will be
!       dumped and cleared, and when nmessg is less than zero the
!       tables will be dumped and not cleared.
!     --output--
!       icount will be the number of times this message has
!       been seen, or zero if the table has overflowed and
!       does not contain this message specifically.
!       when nmessg=0, icount will not be altered.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  19 mar 1980
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer i
  integer i1mach
  integer icount
  integer ii
  integer iunit
  integer kount(10)
  integer, save :: kountx
  integer kunit
  integer level
  integer levtab(10)
  integer lun(5)
  character ( len = 20 ) mes
  character*(*) messg
  character*20 mestab(10)
  integer nerr
  integer, save, dimension ( 10 ) :: nertab
  integer nmessg
  integer nunit
!
  save mestab,levtab,kount
!
!     next two data statements are necessary to provide a blank
!     error table initially
!
  data kount(1),kount(2),kount(3),kount(4),kount(5), &
    kount(6),kount(7),kount(8),kount(9),kount(10) &
    /0,0,0,0,0,0,0,0,0,0/
  data kountx/0/
!
  if (nmessg>0) go to 80
!     dump the table
     if (kount(1)==0) return
!        print to each unit
     call xgetua(lun,nunit)

     do kunit=1,nunit
        iunit = lun(kunit)
        if (iunit==0) iunit = i1mach(4)
!           print table header
        write (iunit,10)
   10       format (32h0          error message summary/ &
            51h message start             nerr     level     count)
!           print body of table
        do i=1,10
           if (kount(i)==0) go to 30
           write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
   15          format (1x,a20,3i10)
        end do
   30       continue
!           print number of other errors
        if (kountx/=0) write (iunit,40) kountx
   40       format (41h0other errors not individually tabulated=,i10)
        write (iunit,50)
   50       format (1x)
     end do

     if (nmessg<0) return
!
!        clear the error tables
!
     do i=1,10
       kount(i) = 0
     end do

     kountx = 0
     return
   80 continue
!     process a message...
!     search for this messg, or else an empty slot for this messg,
!     or else determine that the error table is full.
  mes = messg
  do 90 i=1,10
     ii = i
     if (kount(i)==0) go to 110
     if (mes/=mestab(i)) go to 90
     if (nerr/=nertab(i)) go to 90
     if (level/=levtab(i)) go to 90
     go to 100
   90 continue
!     three possible cases...
!     table is full
     kountx = kountx+1
     icount = 1
     return
!     message found in table
  100    kount(ii) = kount(ii) + 1
     icount = kount(ii)
     return
!     empty slot found for new message
  110    mestab(ii) = mes
     nertab(ii) = nerr
     levtab(ii) = level
     kount(ii)  = 1
     icount = 1
     return
end
subroutine xgetf(kontrl)
!
!*******************************************************************************
!
!! XGETF returns current value of error control flag.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!   abstract
!        xgetf returns the current value of the error control flag
!        in kontrl.  see subroutine xsetf for flag value meanings.
!        (kontrl is an output parameter only.)
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer j4save
  integer kontrl
!
  kontrl = j4save(2,0,.false.)

  return
end
subroutine xgetua(iunita,n)
!
!*******************************************************************************
!
!! XGETUA returns unit number(s) to which error messages are being sent.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xgetua may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        these unit numbers may have been set by a call to xsetun,
!        or a call to xsetua, or may be a default value.
!
!     description of parameters
!      --output--
!        iunit - an array of one to five unit numbers, depending
!                on the value of n.  a value of zero refers to the
!                default unit, as defined by the i1mach machine
!                constant routine.  only iunit(1),...,iunit(n) are
!                defined by xgetua.  the values of iunit(n+1),...,
!                iunit(5) are not defined (for n < 5) or altered
!                in any way by xgetua.
!        n     - the number of units to which copies of the
!                error messages are being sent.  n will be in the
!                range from 1 to 5.
!
!     latest revision ---  19 mar 1980
!     written by ron jones, with slatec common math library subcommittee
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer i
  integer index
  integer iunita(5)
  integer j4save
  integer n
!
  n = j4save(5,0,.false.)

  do i=1,n
     index = i+4
     if (i==1) index = 3
     iunita(i) = j4save(index,0,.false.)
  end do

  return
end
subroutine xgetun(iunit)
!
!*******************************************************************************
!
!! XGETUN returns the (first) output file to which messages are being sent.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xgetun gets the (first) output file to which error messages
!        are being sent.  to find out if more than one file is being
!        used, one must use the xgetua routine.
!
!     description of parameter
!      --output--
!        iunit - the logical unit number of the  (first) unit to
!                which error messages are being sent.
!                a value of zero means that the default file, as
!                defined by the i1mach routine, is being used.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision --- 23 may 1979
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer iunit
  integer j4save
!
  iunit = j4save(3,0,.false.)

  return
end
subroutine xsetf(kontrl)
!
!*******************************************************************************
!
!! XSETF sets the error control flag.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xsetf sets the error control flag value to kontrl.
!        (kontrl is an input parameter only.)
!        the following table shows how each message is treated,
!        depending on the values of kontrl and level.  (see xerror
!        for description of level.)
!
!        if kontrl is zero or negative, no information other than the
!        message itself (including numeric values, if any) will be
!        printed.  if kontrl is positive, introductory messages,
!        trace-backs, etc., will be printed in addition to the message.
!
!              abs(kontrl)
!        level        0              1              2
!        value
!          2        fatal          fatal          fatal
!
!          1     not printed      printed         fatal
!
!          0     not printed      printed        printed
!
!         -1     not printed      printed        printed
!                                  only           only
!                                  once           once
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  19 mar 1980
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer j4save
  integer junk
  integer kontrl
!
  if ((kontrl>=(-2)).and.(kontrl<=2)) go to 10
     call xerrwv('xsetf  -- invalid value of kontrl (i1).',33,1,2, &
       1,kontrl,0,0,0.,0.)
     return
   10 junk = j4save(2,kontrl,.true.)
  return
end
subroutine xsetua(iunita,n)
!
!*******************************************************************************
!
!! XSETUA sets up to 5 unit numbers to which messages are to be sent.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xsetua may be called to declare a list of up to five
!        logical units, each of which is to receive a copy of
!        each error message processed by this package.
!        the purpose of xsetua is to allow simultaneous printing
!        of each error message on, say, a main output file,
!        an interactive terminal, and other files such as graphics
!        communication files.
!
!     description of parameters
!      --input--
!        iunit - an array of up to five unit numbers.
!                normally these numbers should all be different
!                (but duplicates are not prohibited.)
!        n     - the number of unit numbers provided in iunit
!                must have 1 <= n <= 5.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  19 mar 1980
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer i
  integer index
  integer iunita(5)
  integer j4save
  integer junk
  integer n
!
  if ((n>=1).and.(n<=5)) go to 10
     call xerrwv('xsetua -- invalid value of n (i1).',34,1,2,1,n,0,0,0.,0.)
     return
   10 continue

  do i=1,n
     index = i+4
     if (i==1) index = 3
     junk = j4save(index,iunita(i),.true.)
  end do

  junk = j4save(5,n,.true.)

  return
end
subroutine xsetun(iunit)
!
!*******************************************************************************
!
!! XSETUN sets the output file to which error messages are to be sent.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!     abstract
!        xsetun sets the output file to which error messages are to
!        be sent.  only one file will be used.  see xsetua for
!        how to declare more than one file.
!
!     description of parameter
!      --input--
!        iunit - an input parameter giving the logical unit number
!                to which error messages are to be sent.
!
!     written by ron jones, with slatec common math library subcommittee
!     latest revision ---  7 june 1978
!  references  jones r.e., kahaner d.k., "xerror, the slatec error-
!                 handling package", sand82-0800, sandia laboratories,
!                 1982.
!
  integer iunit
  integer j4save
  integer junk
!
  junk = j4save(3,iunit,.true.)
  junk = j4save(5,1,.true.)

  return
end
