!  rkfprb.f90  02 June 2000
!
program rkfprb
!
!*******************************************************************************
!
!! RKFPRB tests the RKF45 ODE integrator.
!
  integer neqn
  parameter ( neqn = 1 )
!
  real abserr
  integer i_step
  integer iflag
  integer iwork(5)
  integer n_step
  real relerr
  real t
  real t_out
  real t_start
  real t_stop
  real work(6*neqn+3)
  real y(neqn)
  real y_exact
!
  external f
!
  abserr = 0.00001
  relerr = 0.00001

  iflag = 1

  t_start = 0.0
  t_stop = 20.0

  n_step = 5

  t_out = 0.0
  y(1) = 1.0

  write ( *, * ) ' '
  write ( *, * ) 'RKFPRB'
  write ( *, * ) '  Use RKF45 to integrate an ODE.'
  write ( *, * ) ' '
  write ( *, * ) '       T          Y         Y_Exact         Error'
  write ( *, * ) ' '
  write ( *, '(4g14.6)' ) t_out, y(1), y_exact(t_out), y(1) - y_exact(t_out)

  do i_step = 1, n_step

    t = ( ( n_step - i_step + 1 ) * t_start &
            + ( i_step - 1 ) * t_stop ) / real ( n_step )
    t_out = ( ( n_step - i_step ) * t_start &
            + ( i_step ) * t_stop ) / real ( n_step )

    call rkf45 ( f, neqn, y, t, t_out, relerr, abserr, iflag, work, iwork )

    write ( *, '(4g14.6)' ) t_out, y(1), y_exact(t_out), y(1) - y_exact(t_out)

  end do

  write ( *, * ) ' '
  write ( *, * ) 'RKFPRB'
  write ( *, * ) '  Normal end of RKF test.'

  stop
end
subroutine f ( t, y, yp )
!
!*******************************************************************************
!
!! F evaluates the derivative for the ODE.
!
  real t
  real y(1)
  real yp(1)
!
  yp(1) = 0.25 * y(1) * ( 1.0 - y(1) / 20.0 )

  return
end
function y_exact ( t )
!
!*******************************************************************************
!
!! Y_EXACT evaluates the exact solution of the ODE.
!
  real t
  real y_exact
!
  y_exact = 20.0 / ( 1.0 + 19.0 * exp ( - 0.25 * t ) )

  return
end
