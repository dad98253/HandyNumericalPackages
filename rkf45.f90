!  rkf45.f90  02 June 2000
!
subroutine rkf45 ( f, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )
!
!*******************************************************************************
!
!! RKF45 carries out the Fehlberg fourth-fifth order Runge-Kutta method.
!
!
!  Author:
!
!    H A Watts and L F Shampine,
!    Sandia Laboratories,
!    Albuquerque, New Mexico.
!
!    RKF45 is primarily designed to solve non-stiff and mildly stiff
!    differential equations when derivative evaluations are inexpensive.
!    RKF45 should generally not be used when the user is demanding
!    high accuracy.
!
!  Abstract:
!
!    RKF45 integrates a system of NEQN first order ordinary differential 
!    equations of the form:
!
!      dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
!
!    where the y(i) are given at t.
!
!    Typically the subroutine is used to integrate from t to tout but it
!    can be used as a one-step integrator to advance the solution a
!    single step in the direction of tout.  on return the parameters in
!    the call list are set for continuing the integration. the user has
!    only to call rkf45 again (and perhaps define a new value for tout).
!    actually, rkf45 is an interfacing routine which calls subroutine
!    rkfs for the solution.  rkfs in turn calls subroutine  fehl which
!    computes an approximate solution over one step.
!
!    RKF45 uses the Runge-Kutta-Fehlberg (4,5)  method described
!    in the reference
!    e.fehlberg , low-order classical runge-kutta formulas with stepsize
!                 control , nasa tr r-315
!
!    the performance of rkf45 is illustrated in the reference
!    l.f.shampine,h.a.watts,s.davenport, solving non-stiff ordinary
!                 differential equations-the state of the art ,
!                 sandia laboratories report sand75-0182 ,
!                 to appear in siam review.
!
!
!  Parameters:
!
!    Input, external F, a subroutine of the form
!      subroutine f(t,y,yp) 
!    to evaluate derivatives yp(i)=dy(i)/dt
!
!    Input, integer NEQN, the number of equations to be integrated.
!
!    Input/output, real Y(NEQN), the solution vector at T.
!
!    Input/output, real T, the independent variable.
!
!    Input, real TOUT, the output point at which solution is desired.
!
!    Input, real RELERR, ABSERR, the relative and absolute error tolerances 
!    for the local error test.  At each step the code requires:
!      abs(local error) <= relerr*abs(y) + abserr
!    for each component of the local error and solution vectors
!
!    Output, integer IFLAG, indicator for status of integration.
!
!    Workspace, real WORK(3+6*NEQN), an array to hold information internal 
!    to RKF45 which is necessary for subsequent calls. 
!
!    Workspace, integer IWORK(5), an array used to hold information internal 
!    to RKF45 which is necessary for subsequent calls. 
!
!
!  first call
!
!    The user must provide storage in his calling program for the arrays
!    in the call list  -      y(neqn) , work(3+6*neqn) , iwork(5)  ,
!    declare f in an external statement, supply subroutine f(t,y,yp) and
!    initialize the following parameters-
!
!      neqn -- number of equations to be integrated.  (neqn >= 1)
!
!      y(*) -- vector of initial conditions
!
!      t -- starting point of integration , must be a variable
!
!      tout -- output point at which solution is desired.
!            t=tout is allowed on the first call only, in which case
!            rkf45 returns with iflag=2 if continuation is possible.
!
!      relerr,abserr -- relative and absolute local error tolerances
!            which must be non-negative. relerr must be a variable while
!            abserr may be a constant. the code should normally not be
!            used with relative error control smaller than about 1.e-8 .
!            to avoid limiting precision difficulties the code requires
!            relerr to be larger than an internally computed relative
!            error parameter which is machine dependent. in particular,
!            pure absolute error is not permitted. if a smaller than
!            allowable value of relerr is attempted, rkf45 increases
!            relerr appropriately and returns control to the user before
!            continuing the integration.
!
!      iflag -- +1,-1  indicator to initialize the code for each new
!            problem. normal input is +1. the user should set iflag=-1
!            only when one-step integrator control is essential. in this
!            case, rkf45 attempts to advance the solution a single step
!            in the direction of tout each time it is called. since this
!            mode of operation results in extra computing overhead, it
!            should be avoided unless needed.
!
!
!  output from rkf45
!
!      y(*) -- solution at t
!      t -- last point reached in integration.
!      iflag = 2 -- integration reached tout. indicates successful retur
!                   and is the normal mode for continuing integration.
!            =-2 -- a single successful step in the direction of tout
!                   has been taken. normal mode for continuing
!                   integration one step at a time.
!            = 3 -- integration was not completed because relative error
!                   tolerance was too small. relerr has been increased
!                   appropriately for continuing.
!            = 4 -- integration was not completed because more than
!                   3000 derivative evaluations were needed. this
!                   is approximately 500 steps.
!            = 5 -- integration was not completed because solution
!                   vanished making a pure relative error test
!                   impossible. must use non-zero abserr to continue.
!                   using the one-step integration mode for one step
!                   is a good way to proceed.
!            = 6 -- integration was not completed because requested
!                   accuracy could not be achieved using smallest
!                   allowable stepsize. user must increase the error
!                   tolerance before continued integration can be
!                   attempted.
!            = 7 -- it is likely that rkf45 is inefficient for solving
!                   this problem. too much output is restricting the
!                   natural stepsize choice. use the one-step integrator
!                   mode.
!            = 8 -- invalid input parameters
!                   this indicator occurs if any of the following is
!                   satisfied -   neqn <= 0
!                                 t=tout  and  iflag /= +1 or -1
!                                 relerr or abserr < 0.
!                                 iflag == 0  or  < -2  or  > 8
!      work(*),iwork(*) -- information which is usually of no interest
!                   to the user but necessary for subsequent calls.
!                   work(1),...,work(neqn) contain the first derivatives
!                   of the solution vector y at t. work(neqn+1) contains
!                   the stepsize h to be attempted on the next step.
!                   iwork(1) contains the derivative evaluation counter.
!
!
!  subsequent calls
!
!    RKF45 returns with all information needed to continue
!    the integration. if the integration reached tout, the user need onl
!    define a new tout and call RKF45 again.  In the one-step integrator
!    mode (iflag=-2) the user must keep in mind that each step taken is
!    in the direction of the current tout.  Upon reaching tout (indicated
!    by changing iflag to 2),the user must then define a new tout and
!    reset iflag to -2 to continue in the one-step integrator mode.
!
!    If the integration was not completed but the user still wants to
!    continue (iflag=3,4 cases), he just calls rkf45 again. with iflag=3
!    the relerr parameter has been adjusted appropriately for continuing
!    the integration. in the case of iflag=4 the function counter will
!    be reset to 0 and another 3000 function evaluations are allowed.
!
!    However,in the case iflag=5, the user must first alter the error
!    criterion to use a positive value of abserr before integration can
!    proceed. if he does not,execution is terminated.
!
!    Also,in the case iflag=6, it is necessary for the user to reset
!    iflag to 2 (or -2 when the one-step integration mode is being used)
!    as well as increasing either abserr,relerr or both before the
!    integration can be continued. if this is not done, execution will
!    be terminated. the occurrence of iflag=6 indicates a trouble spot
!    (solution is changing rapidly,singularity may be present) and it
!    often is inadvisable to continue.
!
!    If iflag=7 is encountered, the user should use the one-step
!    integration mode with the stepsize determined by the code. 
!    If the user insists upon continuing the integration with RKF45, 
!    he must reset iflag to 2 before calling RKF45 again. otherwise,
!    execution will be terminated.
!
!    If iflag=8 is obtained, integration can not be continued unless
!    the invalid input parameters are corrected.
!
!    The arrays work and iwork contain information
!    required for subsequent integration, and should not be altered.
!
  integer neqn
!
  real abserr
  integer iflag
  integer iwork(5)
  integer k1
  integer k1m
  integer k2
  integer k3
  integer k4
  integer k5
  integer k6
  real relerr
  real t
  real tout
  real work(6*neqn+3)
  real y(neqn)
!
  external f
!
!  Compute indices for the splitting of the work array
!
  k1m = neqn + 1
  k1 = k1m + 1
  k2 = k1 + neqn
  k3 = k2 + neqn
  k4 = k3 + neqn
  k5 = k4 + neqn
  k6 = k5 + neqn
!
!  This interfacing routine merely relieves the user of a long
!  calling list via the splitting apart of two working storage arrays.  
!
  call rkfs ( f, neqn, y, t, tout, relerr, abserr, iflag, work(1), work(k1m), &
    work(k1), work(k2), work(k3), work(k4), work(k5), work(k6), work(k6+1), &
    iwork(1), iwork(2), iwork(3), iwork(4), iwork(5) )

  return
end
subroutine rkfs ( f, neqn, y, t, tout, relerr, abserr, iflag, yp, h, &
  f1, f2, f3, f4, f5, savre, savae, &
  nfe, kop, init, jflag, kflag )
!
!*******************************************************************************
!
!! RKFS implements the Fehlberg fourth-fifth order Runge-Kutta method.
!
!
!  Discussion:
!
!    RKFS integrates a system of first order ordinary differential
!    equations as described in the comments for RKF45.
!
!    The arrays yp, f1, f2, f3, f4,and f5 (of dimension at least neqn) and
!    the variables h, savre, savae, nfe, kop, init, jflag and kflag are used
!    internally by the code and appear in the call list to eliminate
!    local retention of variables between calls.  Accordingly, they
!    should not be altered.  Items of possible interest are
!
!      yp - derivative of solution vector at t
!      h  - an appropriate stepsize to be used for the next step
!      nfe- counter on the number of derivative function evaluations
!
!    The expense is controlled by restricting the number
!    of function evaluations to be approximately maxnfe.
!    as set, this corresponds to about 500 steps.
!
!    REMIN is the minimum acceptable value of RELERR.  Attempts
!    to obtain higher accuracy with this subroutine are usually
!    very expensive and often unsuccessful.
!
  real, parameter :: eps = 1.2E-7
  integer, parameter :: maxnfe = 3000
  real, parameter :: remin = 1.0E-12
  real, parameter :: u26 = 26.0 * eps
!
  integer neqn
!
  real a
  real abserr
  real ae
  real dt
  real ee
  real eeoet
  real esttol
  real et
  real f1(neqn)
  real f2(neqn)
  real f3(neqn)
  real f4(neqn)
  real f5(neqn)
  real h
  logical hfaild
  real hmin
  integer iflag
  integer init
  integer jflag
  integer k
  integer kflag
  integer kop
  integer mflag
  integer nfe
  logical output
  real relerr
  real rer
  real s
  real savae
  real savre
  real scale
  real t
  real tol
  real toln
  real tout
  real y(neqn)
  real yp(neqn)
  real ypk
!
  external f
!
!  Check the input parameters.
!
  if ( neqn < 1) then
    iflag = 8
    return
  end if

  if ( relerr < 0.0 ) then
    iflag = 8
    return
  end if

  if ( abserr < 0.0 ) then
    iflag = 8
    return
  end if

  mflag = abs ( iflag )

  if ( abs ( iflag ) < 1 .or. abs ( iflag ) > 8 ) then
    iflag = 8
    return
  end if
!
!  Is this the first call?
!
  if ( mflag == 1 ) then
    go to 50
  end if
!
!  check continuation possibilities
!
  if ( t == tout .and. kflag /= 3 ) then
    iflag = 8
    return
  end if

  if ( mflag /= 2 ) then
    go to 25
  end if
!
!  iflag = +2 or -2
!
  if ( kflag == 3 ) go to 45
  if ( init == 0 ) go to 45
  if ( kflag == 4 ) go to 40

  if ((kflag == 5)  .and.  (abserr == 0.0)) then
    stop
  end if

  if ( kflag == 6 .and. relerr <= savre .and. abserr <= savae ) then
    stop
  end if

  go to 50
!
!  iflag = 3,4,5,6,7 or 8
!
   25 continue

  if ( iflag == 3 ) go to 45
  if ( iflag == 4 ) go to 40
  if ( iflag == 5 .and. abserr > 0.0 ) go to 45
!
!  Integration cannot be continued since user did not respond to
!  the instructions pertaining to iflag=5,6,7 or 8
!
  stop
!
!  Reset function evaluation counter
!
   40 continue

  nfe = 0
  if ( mflag == 2 ) then
    go to 50
  end if
!
!  Reset flag value from previous call
!
   45 continue

  iflag = jflag

  if ( kflag == 3 ) then
    mflag = abs ( iflag )
  end if
!
!  Save input iflag and set continuation flag for subsequent input checking.
!
   50 continue

  jflag = iflag
  kflag = 0
!
!  Save relerr and abserr for checking input on subsequent calls
!
  savre = relerr
  savae = abserr
!
!  Restrict relative error tolerance to be at least as large as
!  2*eps+remin to avoid limiting precision difficulties arising
!  from impossible accuracy requests
!
  rer = 2.0 * epsilon ( rer ) + remin
!
!  The relative error tolerance is too small.
!
  if ( relerr < rer ) then
    relerr = rer
    iflag = 3
    kflag = 3
    return
  end if

  dt = tout - t

  if ( mflag == 1 ) go to 60
  if ( init == 0 ) go to 65
  go to 80
!
!  Initialization:
!    set initialization completion indicator,init
!    set indicator for too many output points,kop
!    evaluate initial derivatives
!    set counter for function evaluations,nfe
!    evaluate initial derivatives
!    set counter for function evaluations,nfe
!    estimate starting stepsize
!
   60 continue

  init = 0
  kop = 0
  a = t
  call f ( a, y, yp )
  nfe = 1

  if ( t == tout ) then
    iflag = 2
    return
  end if

   65 continue

  init = 1
  h = abs ( dt )
  toln = 0.0
  do k = 1, neqn
    tol = relerr * abs ( y(k) ) + abserr
    if ( tol > 0.0 ) then
      toln = tol
      ypk = abs ( yp(k) )
      if ( ypk * h**5 > tol) then
        h = ( tol / ypk )**0.2
      end if
    end if
  end do

  if ( toln <= 0.0 ) then
    h = 0.0
  end if

  h = max ( h, u26 * max ( abs ( t ), abs ( dt ) ) )
  jflag =  sign ( 2, iflag )
!
!  Set stepsize for integration in the direction from t to tout
!
   80 continue

  h = sign ( h, dt )
!
!  Test to see if rkf45 is being severely impacted by too many output points.
!
  if ( abs ( h ) >= 2.0 * abs ( dt ) ) then
    kop = kop + 1
  end if
!
!  Unnecessary frequency of output.
!
  if ( kop == 100 ) then
    kop = 0
    iflag = 7
    return
  end if
!
!  If too close to output point, extrapolate and return.
!
  if ( abs ( dt ) <= u26 * abs ( t ) ) then

    do k = 1, neqn
      y(k) = y(k) + dt * yp(k)
    end do
    a = tout
    call f ( a, y, yp )
    nfe = nfe + 1
    go to 300

  end if
!
!  Initialize output point indicator.
!
  output = .false.
!
!  To avoid premature underflow in the error tolerance function,
!  scale the error tolerances
!
  scale = 2.0 / relerr
  ae = scale * abserr
!
!  Step by step integration.
!
  100 continue

  hfaild = .false.
!
!  Set smallest allowable stepsize.
!
  hmin = u26 * abs ( t )
!
!  Adjust stepsize if necessary to hit the output point.
!  Look ahead two steps to avoid drastic changes in the stepsize and
!  thus lessen the impact of output points on the code.
!
  dt = tout - t
  if ( abs ( dt ) >= 2.0 * abs ( h ) ) go to 200
!
!  The next successful step will complete the integration to the output point.
!
  if ( abs ( dt ) <= abs ( h ) ) then
    output = .true.
    h = dt
    go to 200
  end if

  h = 0.5 * dt
!
!  Core integrator for taking a single step
!
!     the tolerances have been scaled to avoid premature underflow in
!     computing the error tolerance function et.
!     to avoid problems with zero crossings,relative error is measured
!     using the average of the magnitudes of the solution at the
!     beginning and end of a step.
!     the error estimate formula has been grouped to control loss of
!     significance.
!
!     to distinguish the various arguments, h is not permitted
!     to become smaller than 26 units of roundoff in t.
!     practical limits on the change in the stepsize are enforced to
!     smooth the stepsize selection process and to avoid excessive
!     chattering on problems having discontinuities.
!     to prevent unnecessary failures, the code uses 9/10 the stepsize
!     it estimates will succeed.
!
!     after a step failure, the stepsize is not allowed to increase for
!     the next attempted step. this makes the code more efficient on
!     problems having discontinuities and more effective in general
!     since local extrapolation is being used and extra caution seems
!     warranted.
!
!     test number of derivative function evaluations.
!     if okay, try to advance the integration from t to t+h.
!
  200 continue
!
!  Too much work.
!
  if ( nfe > maxnfe ) then
    iflag = 4
    kflag = 4
    return
  end if
!
!  Advance an approximate solution over one step of length H.
!
  call fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, f1 )
  nfe = nfe + 5
!
!  Compute and test allowable tolerances versus local error estimates
!  and remove scaling of tolerances. note that relative error is
!  measured with respect to the average of the magnitudes of the
!  solution at the beginning and end of the step.
!
  eeoet = 0.0
  do k = 1, neqn

    et = abs ( y(k) ) + abs ( f1(k) ) + ae

    if ( et <= 0.0 ) then
      iflag = 5
      return
    end if

    ee = abs ( ( -2090.0 * yp(k) + ( 21970.0 * f3(k) - 15048.0 * f4(k) ) ) + &
      ( 22528.0 * f2(k) - 27360.0 * f5(k) ) )

    eeoet = max ( eeoet, ee / et )

  end do

  esttol = abs ( h ) * eeoet * scale / 752400.0

  if ( esttol <= 1.0) then
    go to 260
  end if
!
!  Unsuccessful step.  Reduce the stepsize, try again.
!  The decrease is limited to a factor of 1/10.
!
  hfaild = .true.
  output = .false.

  if ( esttol < 59049.0 ) then
    s = 0.9 / esttol**0.2
  else
    s = 0.1
  end if

  h = s * h

  if ( abs ( h ) < hmin ) then
    iflag = 6
    kflag = 6
    return
  else
    go to 200
  end if
!
!  Successful step.  Store solution at T+H and evaluate derivatives there.
!
  260 continue

  t = t + h
  do k = 1, neqn
    y(k) = f1(k)
  end do
  a = t
  call f ( a, y, yp )
  nfe = nfe + 1
!
!  Choose next stepsize.  The increase is limited to a factor of 5.
!  If step failure has just occurred, next stepsize is not allowed to increase
!
  if ( esttol > 0.0001889568 ) then
    s = 0.9 / esttol**0.2
  else
    s = 5.0
  end if

  if ( hfaild ) then
    s = min ( s, 1.0 )
  end if

  h = sign ( max ( s * abs ( h ), hmin ), h )
!
!  end of core integrator
!
!  should we take another step
!
  if ( output ) go to 300
  if ( iflag > 0 ) go to 100
!
!  integration successfully completed
!
!  one-step mode
!
  iflag = - 2
  return
!
!  interval mode
!
  300 continue

  t = tout
  iflag = 2

  return
end
subroutine fehl ( f, neqn, y, t, h, yp, f1, f2, f3, f4, f5, s )
!
!*******************************************************************************
!
!! FEHL takes a single Fehlberg fourth-fifth order Runge-Kutta step.
!
!
!  Discussion:
!
!    FEHL integrates a system of NEQN first order ordinary differential 
!    equations of the form
!      dY(i)/dT = F(T,Y(1),---,Y(NEQN))
!    where the initial values Y and the initial derivatives
!    YP are specified at the starting point T.  FEHL advances
!    the solution over the fixed step H and returns
!    the fifth order (sixth order accurate locally) solution
!    approximation at T+H in array S.
!
!    F1,F2,f3,f4,f5 are arrays of dimension NEQN which are needed
!    for internal storage.
!
!    The formulas have been grouped to control loss of significance.
!    FEHL should be called with an H not smaller than 13 units of
!    roundoff in T so that the various independent arguments can be
!    distinguished.
!
  integer neqn
!
  real ch
  real f1(neqn)
  real f2(neqn)
  real f3(neqn)
  real f4(neqn)
  real f5(neqn)
  real h
  integer k
  real s(neqn)
  real t
  real y(neqn)
  real yp(neqn)
!
  ch = h / 4.0
  do k = 1, neqn
    f5(k) = y(k) + ch * yp(k)
  end do
  call f ( t + ch, f5, f1 )

  ch = 3.0 * h / 32.0
  do k = 1, neqn
    f5(k) = y(k) + ch * ( yp(k) + 3.0 * f1(k) )
  end do
  call f ( t + 3.0 * h / 8.0, f5, f2 )

  ch = h / 2197.0
  do k = 1, neqn
    f5(k) = y(k) + ch * ( 1932.0 * yp(k) + ( 7296.0 * f2(k) - 7200.0 * f1(k) ) )
  end do
  call f ( t + 12.0 * h / 13.0, f5, f3 )

  ch = h / 4104.0
  do k = 1, neqn
    f5(k) = y(k) + ch * ( ( 8341.0 * yp(k) - 845.0 * f3(k) ) &
      + ( 29440.0 * f2(k) - 32832.0 * f1(k) ) )
  end do
  call f ( t + h, f5, f4 )

  ch = h / 20520.0
  do k = 1, neqn
    f1(k) = y(k) + ch * ( ( -6080.0 * yp(k) + ( 9295.0 * f3(k) &
      - 5643.0 * f4(k) ) ) + ( 41040.0 * f1(k) - 28352.0 * f2(k) ) )
  end do
  call f ( t + h / 2.0, f1, f5 )
!
!  Compute approximate solution at t+h
!
  ch = h / 7618050.0
  do k = 1, neqn
    s(k) = y(k) + ch * ( ( 902880.0 * yp(k) + ( 3855735.0 * f3(k) &
      - 1371249.0 * f4(k) ) ) + ( 3953664.0 * f2(k) + 277020.0 * f5(k) ) )
  end do

  return
end
