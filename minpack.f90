!  minpack.f90  17 June 2000
!
subroutine chkder ( m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err )
!
!*******************************************************************************
!
!! CHKDER checks the gradients of M functions of N variables.
!
!
!     CHKDER checks the gradients of m nonlinear functions
!     in n variables, evaluated at a point x, for consistency with
!     the functions themselves. the user must call chkder twice,
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
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
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
!         then if err(i) is 1.0 the i-th gradient is correct,
!         while if err(i) is 0.0 the i-th gradient is incorrect.
!         for values of err between 0.0 and 1.0, the categorization
!         is less certain. in general, a value of err(i) greater
!         than 0.5 indicates that the i-th gradient is probably
!         correct, while a value of err(i) less than 0.5 indicates
!         that the i-th gradient is probably incorrect.
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
  real fjac(ldfjac,n)
  real fvec(m)
  real fvecp(m)
  integer i
  integer j
  integer mode
  real temp
  real x(n)
  real xp(n)
!
  epsmch = epsilon ( 1.0 )
  eps = sqrt ( epsmch )
!
!  MODE = 1.
!
  if ( mode == 1 ) then

     do j = 1, n
       temp = eps * abs ( x(j) )
       if (temp == 0.0 ) temp = eps
       xp(j) = x(j) + temp
     end do
!
!  MODE = 2.
!
  else if ( mode == 2 ) then

     epsf = 100.0 * epsmch
     epslog = log10 ( eps )
 
     err = 0.0
 
     do j = 1, n
       temp = abs ( x(j) )
       if ( temp == 0.0 ) temp = 1.0
       do i = 1, m
         err(i) = err(i) + temp * fjac(i,j)
       end do
     end do

     do i = 1, m

       temp = 1.0

       if ( fvec(i) /= 0.0 .and. fvecp(i) /= 0.0 .and. &
         abs(fvecp(i)-fvec(i)) >= epsf*abs(fvec(i)) ) then
         temp = eps * abs((fvecp(i)-fvec(i))/eps-err(i)) &
           /(abs(fvec(i)) + abs(fvecp(i)))
       end if

       err(i) = 1.0

       if ( temp > epsmch .and. temp < eps ) then
         err(i) = ( log10 ( temp ) - epslog ) / epslog
       end if

       if ( temp >= eps ) then
         err(i) = 0.0
       end if

     end do

  end if

  return
end
subroutine dogleg ( n, r, lr, diag, qtb, delta, x )
!
!*******************************************************************************
!
!! DOGLEG finds the minimizing combination of Gauss-Newton and gradient steps.
!
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
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
  integer lr
  integer n
!
  real alpha
  real bnorm
  real delta
  real diag(n)
  real enorm
  real epsmch
  real gnorm
  integer i
  integer j
  integer jj
  integer k
  integer l
  real qnorm
  real qtb(n)
  real r(lr)
  real sgnorm
  real sum
  real temp
  real wa1(n)
  real wa2(n)
  real x(n)
!
  epsmch = epsilon ( 1.0 )
!
!  Calculate the Gauss-Newton direction.
!
  jj = ( n * ( n + 1 ) ) / 2 + 1

  do k = 1, n

     j = n - k + 1
     jj = jj - k
     l = jj + 1
     sum = 0.0

     do i = j+1, n
       sum = sum + r(l)*x(i)
       l = l + 1
     end do

     temp = r(jj)

     if ( temp == 0.0 ) then

       l = j
       do i = 1, j
         temp = max ( temp, abs(r(l)) )
         l = l + n - i
       end do

       if ( temp == 0.0 ) then
         temp = epsmch
       else
         temp = epsmch * temp
       end if

     end if

     x(j) = ( qtb(j) - sum ) / temp

  end do
!
!  Test whether the Gauss-Newton direction is acceptable.
!
  wa1(1:n) = 0.0
  wa2(1:n) = diag(1:n) * x(1:n)
  qnorm = enorm ( n, wa2 )

  if ( qnorm <= delta ) then
    return
  end if
!
!  The Gauss-Newton direction is not acceptable.
!  Calculate the scaled gradient direction.
!
  l = 1
  do j = 1, n
     temp = qtb(j)
     do i = j, n
        wa1(i) = wa1(i) + r(l)*temp
        l = l + 1
     end do
     wa1(j) = wa1(j) / diag(j)
  end do
!
!  Calculate the norm of the scaled gradient.
!  Test for the special case in which the scaled gradient is zero.
!
  gnorm = enorm(n,wa1)
  sgnorm = 0.0
  alpha = delta / qnorm

  if ( gnorm /= 0.0 ) then
!
!  Calculate the point along the scaled gradient which minimizes the quadratic.
!
    wa1(1:n) = ( wa1(1:n) / gnorm ) / diag(1:n)

    l = 1
    do j = 1, n
      sum = 0.0
      do i = j, n
        sum = sum + r(l)*wa1(i)
        l = l + 1
      end do
      wa2(j) = sum
    end do

    temp = enorm ( n, wa2 )
    sgnorm = (gnorm/temp)/temp
!
!  Test whether the scaled gradient direction is acceptable.
!
    alpha = 0.0
!
!  The scaled gradient direction is not acceptable.
!  Calculate the point along the dogleg at which the quadratic is minimized.
!
    if ( sgnorm < delta ) then

      bnorm = enorm(n,qtb)
      temp = (bnorm/gnorm) * (bnorm/qnorm) * (sgnorm/delta)
      temp = temp - (delta/qnorm)*(sgnorm/delta)**2 &
        + sqrt((temp-(delta/qnorm))**2 &
        +(1.0-(delta/qnorm)**2)*(1.0-(sgnorm/delta)**2))

      alpha = ((delta/qnorm)*(1.0 - (sgnorm/delta)**2))/temp

    end if

  end if
!
!  Form appropriate convex combination of the Gauss-Newton
!  direction and the scaled gradient direction.
!
  temp = ( 1.0 - alpha ) * min ( sgnorm, delta )

  x(1:n) = temp * wa1(1:n) + alpha * x(1:n)

  return
end
function enorm ( n, x )
!
!*******************************************************************************
!
!! ENORM computes the Euclidean norm of a vector.
!
!
!     the euclidean norm is computed by accumulating the sum of
!     squares in three different sums. the sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. non-destructive underflows are permitted. underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     the definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant. the main
!     restrictions on these constants are that rdwarf**2 not
!     underflow and rgiant**2 not overflow. the constants
!     given here are suitable for every known computer.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
  integer n
!
  real agiant
  real enorm
  integer i
  real rdwarf
  real rgiant
  real s1
  real s2
  real s3
  real x(n)
  real xabs
  real x1max
  real x3max
!
  rdwarf = sqrt ( tiny ( 1.0 ) )
  rgiant = sqrt ( huge ( 1.0 ) )

  s1 = 0.0
  s2 = 0.0
  s3 = 0.0
  x1max = 0.0
  x3max = 0.0
  agiant = rgiant / real ( n )

  do i = 1, n

    xabs = abs ( x(i) )

    if ( xabs <= rdwarf ) then

      if ( xabs > x3max ) then
        s3 = 1.0 + s3 * ( x3max / xabs )**2
        x3max = xabs
      else if ( xabs /= 0.0 ) then
        s3 = s3 + ( xabs / x3max )**2
      end if

    else if ( xabs >= agiant ) then

      if ( xabs > x1max ) then
        s1 = 1.0 + s1 * ( x1max / xabs )**2
        x1max = xabs
      else
        s1 = s1 + ( xabs / x1max )**2
      end if

    else

      s2 = s2 + xabs**2
   
    end if

  end do
!
!  calculation of norm.
!
  if ( s1 /= 0.0 ) then

    enorm = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )

  else if ( s2 /= 0.0 ) then

    if ( s2 >= x3max ) then
      enorm = sqrt ( s2 * ( 1.0 + ( x3max / s2 ) * ( x3max * s3 ) ) )
    else
      enorm = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
    end if

  else

    enorm = x3max * sqrt ( s3 )

  end if

  return
end
subroutine fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn )
!
!*******************************************************************************
!
!! FDJAC1 estimates an N by N jacobian matrix using forward differences.
!
!
!     this subroutine computes a forward-difference approximation
!     to the n by n jacobian matrix associated with a specified
!     problem of n functions in n variables. if the jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an input array of length n.
!
!       fvec is an input array of length n which must contain the
!         functions evaluated at x.
!
!       fjac is an output n by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac1. see description of fcn.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
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
  real temp
  real wa1(n)
  real wa2(n)
  real x(n)
!
  external fcn
!
  epsmch = epsilon ( 1.0 )

  eps = sqrt ( max ( epsfcn, epsmch ) )
  msum = ml + mu + 1
!
!  Computation of dense approximate jacobian.
!
  if ( msum >= n ) then

     do j = 1, n

        temp = x(j)
        h = eps * abs ( temp )
        if (h == 0.0) h = eps
        x(j) = temp + h
        call fcn(n,x,wa1,iflag)

        if (iflag < 0) then
          exit
        end if

        x(j) = temp
        fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h

     end do

  else
!
!  computation of banded approximate jacobian.
!
     do k = 1, msum

        do j = k, n, msum
           wa2(j) = x(j)
           h = eps*abs(wa2(j))
           if (h == 0.0) h = eps
           x(j) = wa2(j) + h
        end do

        call fcn(n,x,wa1,iflag)

        if (iflag < 0) then
          exit
        end if

        do j = k, n, msum

           x(j) = wa2(j)
           h = eps*abs(wa2(j))
           if (h == 0.0) h = eps

           fjac(1:n,j) = 0.0

           do i = 1, n
              if (i >= j - mu .and. i <= j + ml) then
                fjac(i,j) = (wa1(i) - fvec(i))/h
              end if
           end do

        end do

     end do

  end if

  return
end
subroutine fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )
!
!*******************************************************************************
!
!! FDJAC2 estimates an M by N jacobian matrix using forward differences.
!
!
!     this subroutine computes a forward-difference approximation
!     to the m by n jacobian matrix associated with a specified
!     problem of m functions in n variables.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         real x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!       end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac2.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an input array of length n.
!
!       fvec is an input array of length m which must contain the
!         functions evaluated at x.
!
!       fjac is an output m by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac2. see description of fcn.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
  integer ldfjac
  integer m
  integer n
!
  real eps
  real epsfcn
  real epsmch
  real fjac(ldfjac,n)
  real fvec(m)
  real h
  integer i
  integer iflag
  integer j
  real temp
  real wa(m)
  real x(n)
!
  external fcn
!
  epsmch = epsilon ( 1.0 )
!
  eps = sqrt ( max ( epsfcn, epsmch ) )

  do j = 1, n

     temp = x(j)
     h = eps*abs(temp)
     if (h == 0.0) h = eps
     x(j) = temp + h
     call fcn(m,n,x,wa,iflag)

     if ( iflag < 0 ) then
       exit
     end if

     x(j) = temp
     fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

  end do

  return
end
subroutine hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
  factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf )
!
!*******************************************************************************
!
!! HYBRD seeks a zero of N nonlinear equations in N variables using Powell's method.
!
!
!     the purpose of hybrd is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least maxfev
!         by the end of an iteration.
!
!       ml is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         ml to at least n - 1.
!
!       mu is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         jacobian matrix. if the jacobian is not banded, set
!         mu to at least n - 1.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   relative error between two consecutive iterates
!                    is at most xtol.
!
!         info = 2   number of calls to fcn has reached or exceeded
!                    maxfev.
!
!         info = 3   xtol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    five jacobian evaluations.
!
!         info = 5   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    ten iterations.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       r is an output array of length lr which contains the
!         upper triangular matrix produced by the qr factorization
!         of the final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains
!         the vector (q transpose)*fvec.
!
  integer ldfjac
  integer lr
  integer n
!
  real actred
  real delta
  real diag(n)
  real enorm
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
  integer iter
  integer iwa(1)
  integer j
  logical jeval
  integer l
  integer maxfev
  integer ml
  integer mode
  integer msum
  integer mu
  integer ncfail
  integer nslow1
  integer nslow2
  integer ncsuc
  integer nfev
  integer nprint
  real pnorm
  real prered
  real qtf(n)
  real r(lr)
  real ratio
  logical sing
  real sum
  real temp
  real wa1(n)
  real wa2(n)
  real wa3(n)
  real wa4(n)
  real x(n)
  real xnorm
  real xtol
!
  external fcn
!
  epsmch = epsilon ( 1.0 )

  info = 0
  iflag = 0
  nfev = 0
!
!  check the input parameters for errors.
!
  if ( n <= 0 ) then
    return
  else if ( xtol < 0.0 ) then
    return
  else if ( maxfev <= 0 ) then
    return
  else if ( ml < 0 ) then
    return
  else if ( mu < 0 ) then
    return
  else if ( factor <= 0.0 ) then
    return
  else if ( ldfjac < n ) then
    return
  else if ( lr < ( n * (n + 1) ) / 2 ) then
    return
  end if

  if ( mode == 2 ) then

    do j = 1, n
      if (diag(j) <= 0.0) go to 300
    end do

  end if
!
!  evaluate the function at the starting point
!  and calculate its norm.
!
  iflag = 1
  call fcn(n,x,fvec,iflag)
  nfev = 1
  if (iflag < 0) go to 300
  fnorm = enorm(n,fvec)
!
!  determine the number of calls to fcn needed to compute the jacobian matrix.
!
  msum = min ( ml+mu+1,n)
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
     iflag = 2
     call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
     nfev = nfev + msum
     if (iflag < 0) go to 300
!
!  compute the qr factorization of the jacobian.
!
     call qrfac ( n, n, fjac, ldfjac, .false., iwa, 1, wa1, wa2 )
!
!  on the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
    if ( iter == 1 ) then

      if ( mode /= 2 ) then

        diag(1:n) = wa2(1:n)
        do j = 1, n
          if (wa2(j) == 0.0) diag(j) = 1.0
        end do

      end if
!
!  on the first iteration, calculate the norm of the scaled x
!  and initialize the step bound delta.
!
      wa3(1:n) = diag(1:n) * x(1:n)
      xnorm = enorm(n,wa3)
      delta = factor*xnorm
      if (delta == 0.0) delta = factor

    end if
!
!  form (q transpose)*fvec and store in qtf.
!
     qtf(1:n) = fvec(1:n)

     do j = 1, n

       if ( fjac(j,j) /= 0.0 ) then

         sum = 0.0
         do i = j, n
           sum = sum + fjac(i,j) * qtf(i)
         end do

         temp = - sum / fjac(j,j)
         qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp

       end if

     end do
!
!  copy the triangular factor of the qr factorization into r.
!
     sing = .false.

     do j = 1, n
        l = j
        do i = 1, j-1
           r(l) = fjac(i,j)
           l = l + n - i
        end do
        r(l) = wa1(j)
        if (wa1(j) == 0.0) sing = .true.
     end do
!
!  Accumulate the orthogonal factor in FJAC.
!
     call qform ( n, n, fjac, ldfjac )
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  Beginning of the inner loop.
!
  180    continue
!
!  if requested, call fcn to enable printing of iterates.
!
        if (nprint > 0) then
          iflag = 0
          if (mod(iter-1,nprint) == 0) then
            call fcn(n,x,fvec,iflag)
          end if
          if (iflag < 0) go to 300
        end if
!
!  Determine the direction P.
!
        call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
!
!  store the direction p and x + p. calculate the norm of p.
!
        wa1(1:n) = - wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm ( n, wa3 )
!
!  on the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  evaluate the function at x + p and calculate its norm.
!
        iflag = 1
        call fcn(n,wa2,wa4,iflag)
        nfev = nfev + 1
        if (iflag < 0) go to 300
        fnorm1 = enorm(n,wa4)
!
!  compute the scaled actual reduction.
!
        actred = -1.0
        if (fnorm1 < fnorm) actred = 1.0 - (fnorm1/fnorm)**2
!
!  compute the scaled predicted reduction.
!
        l = 1
        do i = 1, n
           sum = 0.0
           do j = i, n
              sum = sum + r(l)*wa1(j)
              l = l + 1
           end do
           wa3(i) = qtf(i) + sum
        end do

        temp = enorm(n,wa3)
        prered = 0.0
        if (temp < fnorm) prered = 1.0 - (temp/fnorm)**2
!
!  compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0
        if (prered > 0.0) ratio = actred/prered
!
!  update the step bound.
!
        if ( ratio < 0.1 ) then
           ncsuc = 0
           ncfail = ncfail + 1
           delta = 0.5 * delta
        else
           ncfail = 0
           ncsuc = ncsuc + 1
           if (ratio >= 0.5 .or. ncsuc > 1) then
             delta = max ( delta,pnorm/ 0.5 )
           end if
           if (abs(ratio-1.0) <= 0.1 ) delta = pnorm/ 0.5
        end if
!
!  test for successful iteration.
!

!
!  successful iteration. update x, fvec, and their norms.
!
        if ( ratio >= 0.0001 ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:n) = wa4(1:n)
          xnorm = enorm(n,wa2)
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  determine the progress of the iteration.
!
        nslow1 = nslow1 + 1
        if (actred >= 0.001 ) nslow1 = 0
        if (jeval) nslow2 = nslow2 + 1
        if (actred >= 0.1 ) nslow2 = 0
!
!  test for convergence.
!
        if (delta <= xtol*xnorm .or. fnorm == 0.0) info = 1
        if (info /= 0) go to 300
!
!  tests for termination and stringent tolerances.
!
        if (nfev >= maxfev) info = 2
        if (0.1 *max ( 0.1 *delta,pnorm) <= epsmch*xnorm) info = 3
        if (nslow2 == 5) info = 4
        if (nslow1 == 10) info = 5
        if (info /= 0) go to 300
!
!  criterion for recalculating jacobian approximation
!  by forward differences.
!
        if (ncfail == 2) go to 290
!
!  calculate the rank one modification to the jacobian
!  and update qtf if necessary.
!
        do j = 1, n
           sum = 0.0
           do i = 1, n
              sum = sum + fjac(i,j)*wa4(i)
           end do
           wa2(j) = (sum - wa3(j))/pnorm
           wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
           if (ratio >= 0.0001 ) qtf(j) = sum
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
!      end of the outer loop.
!
     go to 30
  300 continue
!
!  termination, either normal or user imposed.
!
  if (iflag < 0) info = iflag
  iflag = 0

  if ( nprint > 0 ) then
    call fcn(n,x,fvec,iflag)
  end if

  return
end
subroutine hybrd1 ( fcn, n, x, fvec, tol, info )
!
!*******************************************************************************
!
!! HYBRD1 seeks a zero of N nonlinear equations in N variables using Powell's method.
!
!
!     the purpose of hybrd1 is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. this is done by using the
!     more general nonlinear equation solver hybrd. the user
!     must provide a subroutine which calculates the functions.
!     the jacobian is then calculated by a forward-difference
!     approximation.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrd1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates that the relative error
!         between x and the solution is at most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   algorithm estimates that the relative error
!                    between x and the solution is at most tol.
!
!         info = 2   number of calls to fcn has reached or exceeded
!                    200*(n+1).
!
!         info = 3   tol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress.
!
  integer lwa
  integer n
!
  real diag(n)
  real epsfcn
  real factor
  real fjac(n,n)
  real fvec(n)
  integer info
  integer j
  integer ldfjac
  integer lr
  integer maxfev
  integer ml
  integer mode
  integer mu
  integer nfev
  integer nprint
  real qtf(n)
  real r((n*(n+1))/2)
  real tol
  real x(n)
  real xtol
!
  external fcn
!
  info = 0

  if ( n <= 0 ) then
    return
  else if ( tol < 0.0 ) then
    return
  end if

  maxfev = 200 * ( n + 1 )
  xtol = tol
  ml = n - 1
  mu = n - 1
  epsfcn = 0.0
  mode = 2
  diag(1:n) = 1.0
  nprint = 0
  lr = ( n * ( n + 1 ) ) / 2
  factor = 100.0
  ldfjac = n

  call hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
    factor, nprint, info, nfev, fjac, ldfjac, r, lr, qtf )

  if ( info == 5 ) then
    info = 4
  end if

  return
end
subroutine hybrj ( fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
  factor, nprint, info, nfev, njev, r, lr, qtf )
!
!*******************************************************************************
!
!! HYBRJ seeks a zero of N nonlinear equations in N variables by Powell's method.
!
!
!     the purpose of hybrj is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         real x(n),fvec(n),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ---------
!         return
!       end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrj.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. fvec and fjac should not be altered.
!         if nprint is not positive, no special calls of fcn
!         with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   relative error between two consecutive iterates
!                    is at most xtol.
!
!         info = 2   number of calls to fcn with iflag = 1 has
!                    reached maxfev.
!
!         info = 3   xtol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    five jacobian evaluations.
!
!         info = 5   iteration is not making good progress, as
!                    measured by the improvement from the last
!                    ten iterations.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       r is an output array of length lr which contains the
!         upper triangular matrix produced by the qr factorization
!         of the final approximate jacobian, stored rowwise.
!
!       lr is a positive integer input variable not less than
!         (n*(n+1))/2.
!
!       qtf is an output array of length n which contains
!         the vector (q transpose)*fvec.
!
  integer ldfjac
  integer lr
  integer n
!
  real actred
  real delta
  real diag(n)
  real enorm
  real epsmch
  real factor
  real fjac(ldfjac,n)
  real fnorm
  real fnorm1
  real fvec(n)
  integer i
  integer iflag
  integer info
  integer iter
  integer iwa(1)
  integer j
  logical jeval
  integer l
  integer maxfev
  integer mode
  integer ncfail
  integer nslow1
  integer nslow2
  integer ncsuc
  integer nfev
  integer njev
  integer nprint
  real pnorm
  real prered
  real qtf(n)
  real r(lr)
  real ratio
  logical sing
  real sum
  real temp
  real wa1(n)
  real wa2(n)
  real wa3(n)
  real wa4(n)
  real x(n)
  real xnorm
  real xtol
!
  external fcn
!
  epsmch = epsilon ( 1.0 )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  check the input parameters for errors.
!
  if (n <= 0 .or. ldfjac < n .or. xtol < 0.0 &
    .or. maxfev <= 0 .or. factor <= 0.0 &
    .or. lr < (n*(n + 1))/2) go to 300

  if ( mode == 2 ) then
    do j = 1, n
      if (diag(j) <= 0.0) go to 300
    end do
  end if

   20 continue
!
!  evaluate the function at the starting point
!  and calculate its norm.
!
  iflag = 1
  call fcn(n,x,fvec,fjac,ldfjac,iflag)
  nfev = 1
  if (iflag < 0) go to 300
  fnorm = enorm(n,fvec)
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
     iflag = 2
     call fcn(n,x,fvec,fjac,ldfjac,iflag)
     njev = njev + 1
     if (iflag < 0) go to 300
!
!  compute the qr factorization of the jacobian.
!
     call qrfac ( n, n, fjac, ldfjac, .false., iwa, 1, wa1, wa2 )
!
!  on the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if ( iter == 1) then

       if (mode /= 2) then
         diag(1:n) = wa2(1:n)
         do j = 1, n
           if (wa2(j) == 0.0) diag(j) = 1.0
         end do
       end if
!
!  on the first iteration, calculate the norm of the scaled x
!  and initialize the step bound delta.
!
       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm(n,wa3)
       delta = factor*xnorm
       if (delta == 0.0) delta = factor
     end if
!
!  form (q transpose)*fvec and store in qtf.
!
     qtf(1:n) = fvec(1:n)

     do j = 1, n
       if ( fjac(j,j) /= 0.0 ) then
         sum = 0.0
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
        do i = 1, j-1
           r(l) = fjac(i,j)
           l = l + n - i
        end do
        r(l) = wa1(j)
        if (wa1(j) == 0.0) sing = .true.
     end do
!
!  Accumulate the orthogonal factor in FJAC.
!
     call qform ( n, n, fjac, ldfjac )
!
!  Rescale if necessary.
!
     if (mode /= 2) then
       do j = 1, n
         diag(j) = max ( diag(j),wa2(j))
       end do
     end if
!
!  beginning of the inner loop.
!
  180    continue
!
!  if requested, call fcn to enable printing of iterates.
!
        if (nprint > 0) then
          iflag = 0
          if (mod(iter-1,nprint) == 0) then
            call fcn(n,x,fvec,fjac,ldfjac,iflag)
          end if
          if (iflag < 0) go to 300
        end if
!
!  Determine the direction P.
!
        call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
!
!  store the direction p and x + p. calculate the norm of p.
!
        wa1(1:n) = - wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)
        pnorm = enorm(n,wa3)
!
!  on the first iteration, adjust the initial step bound.
!
        if (iter == 1) delta = min ( delta,pnorm)
!
!  evaluate the function at x + p and calculate its norm.
!
        iflag = 1
        call fcn(n,wa2,wa4,fjac,ldfjac,iflag)
        nfev = nfev + 1
        if (iflag < 0) go to 300
        fnorm1 = enorm(n,wa4)
!
!  compute the scaled actual reduction.
!
        actred = -1.0
        if (fnorm1 < fnorm) actred = 1.0 - (fnorm1/fnorm)**2
!
!  compute the scaled predicted reduction.
!
        l = 1
        do i = 1, n
           sum = 0.0
           do j = i, n
              sum = sum + r(l)*wa1(j)
              l = l + 1
           end do
           wa3(i) = qtf(i) + sum
        end do

        temp = enorm(n,wa3)
        prered = 0.0
        if (temp < fnorm) prered = 1.0 - (temp/fnorm)**2
!
!  compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0
        if (prered > 0.0) ratio = actred/prered
!
!  update the step bound.
!
        if ( ratio < 0.1 ) then
           ncsuc = 0
           ncfail = ncfail + 1
           delta = 0.5 * delta
        else
           ncfail = 0
           ncsuc = ncsuc + 1

           if (ratio >= 0.5 .or. ncsuc > 1) then
             delta = max ( delta,pnorm/0.5)
           end if

           if (abs(ratio-1.0) <= 0.1 ) delta = pnorm/0.5
        end if
!
!  test for successful iteration.
!

!
!  successful iteration. update x, fvec, and their norms.
!
        if (ratio >= 0.0001 ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:n) = wa4(1:n)
          xnorm = enorm(n,wa2)
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  determine the progress of the iteration.
!
        nslow1 = nslow1 + 1
        if (actred >= 0.001 ) nslow1 = 0
        if (jeval) nslow2 = nslow2 + 1
        if (actred >= 0.1) nslow2 = 0
!
!  test for convergence.
!
        if (delta <= xtol*xnorm .or. fnorm == 0.0) info = 1
        if (info /= 0) go to 300
!
!  tests for termination and stringent tolerances.
!
        if (nfev >= maxfev) info = 2
        if (0.1 *max ( 0.1 *delta,pnorm) <= epsmch*xnorm) info = 3
        if (nslow2 == 5) info = 4
        if (nslow1 == 10) info = 5
        if (info /= 0) go to 300
!
!  criterion for recalculating jacobian.
!
        if (ncfail == 2) go to 290
!
!  calculate the rank one modification to the jacobian
!  and update qtf if necessary.
!
        do j = 1, n
           sum = 0.0
           do i = 1, n
              sum = sum + fjac(i,j)*wa4(i)
           end do
           wa2(j) = (sum - wa3(j))/pnorm
           wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
           if (ratio >= 0.0001 ) qtf(j) = sum
        end do
!
!  compute the qr factorization of the updated jacobian.
!
        call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
        call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
        call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
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

  if ( nprint > 0 ) then
    call fcn(n,x,fvec,fjac,ldfjac,iflag)
  end if

  return
end
subroutine hybrj1 ( fcn, n, x, fvec, fjac, ldfjac, tol, info )
!
!*******************************************************************************
!
!! HYBRJ1 seeks a zero of N nonlinear equations in N variables by Powell's method.
!
!
!     the purpose of hybrj1 is to find a zero of a system of
!     n nonlinear functions in n variables by a modification
!     of the powell hybrid method. this is done by using the
!     more general nonlinear equation solver hybrj. the user
!     must provide a subroutine which calculates the functions
!     and the jacobian.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!         integer n,ldfjac,iflag
!         real x(n),fvec(n),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ---------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of hybrj1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length n which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array which contains the
!         orthogonal matrix q produced by the qr factorization
!         of the final approximate jacobian.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates that the relative error
!         between x and the solution is at most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0   improper input parameters.
!
!         info = 1   algorithm estimates that the relative error
!                    between x and the solution is at most tol.
!
!         info = 2   number of calls to fcn with iflag = 1 has
!                    reached 100*(n+1).
!
!         info = 3   tol is too small. no further improvement in
!                    the approximate solution x is possible.
!
!         info = 4   iteration is not making good progress.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         (n*(n+13))/2.
!
  integer ldfjac
  integer n
!
  real diag(n)
  real factor
  real fjac(ldfjac,n)
  real fvec(n)
  integer info
  integer j
  integer lr
  integer maxfev
  integer mode
  integer nfev
  integer njev
  integer nprint
  real qtf(n)
  real r((n*(n+1))/2)
  real tol
  real x(n)
  real xtol
!
  external fcn
!
  info = 0

  if ( n <= 0 ) then
    return
  else if ( ldfjac < n ) then
    return 
  else if ( tol < 0.0 ) then
    return
  end if

  maxfev = 100 * ( n + 1 )
  xtol = tol
  mode = 2
  diag(1:n) = 1.0
  factor = 100.0
  nprint = 0
  lr = ( n * ( n + 1 ) ) / 2

  call hybrj ( fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, &
    factor, nprint, info, nfev, njev, r, lr, qtf )

  if ( info == 5 ) then
    info = 4
  end if

  return
end
subroutine lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
  diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )
!
!*******************************************************************************
!
!! LMDER minimizes M functions in N variables using the Levenberg-Marquardt method.
!
!
!     the purpose of LMDER is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the Levenberg-Marquardt algorithm. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         real x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of LMDER.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x, fvec, and fjac
!         available for printing. fvec and fjac should not be
!         altered. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
  integer ldfjac
  integer m
  integer n
!
  real actred
  real delta
  real diag(n)
  real dirder
  real enorm
  real epsmch
  real factor
  real fjac(ldfjac,n)
  real fnorm
  real fnorm1
  real ftol
  real fvec(m)
  real gnorm
  real gtol
  integer i
  integer iflag
  integer info
  integer ipvt(n)
  integer iter
  integer j
  integer l
  integer maxfev
  integer mode
  integer nfev
  integer njev
  integer nprint
  real par
  real pnorm
  real prered
  real qtf(n)
  real ratio
  real sum
  real temp
  real temp1
  real temp2
  real wa1(n)
  real wa2(n)
  real wa3(n)
  real wa4(m)
  real xnorm
  real x(n)
  real xtol
!
  external fcn
!
  epsmch = epsilon ( 1.0 )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 .or. m < n .or. ldfjac < m &
    .or. ftol < 0.0 .or. xtol < 0.0 .or. gtol < 0.0 &
     .or. maxfev <= 0 .or. factor <= 0.0) go to 300

  if ( mode == 2 ) then
    do j = 1, n
      if (diag(j) <= 0.0) go to 300
    end do
  end if
!
!  evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
  nfev = 1
  if (iflag < 0) go to 300
  fnorm = enorm(m,fvec)
!
!  initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0
  iter = 1
!
!  beginning of the outer loop.
!
   30 continue
!
!  calculate the jacobian matrix.
!
     iflag = 2
     call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
     njev = njev + 1
     if (iflag < 0) go to 300
!
!  if requested, call fcn to enable printing of iterates.
!
     if ( nprint > 0 ) then
       iflag = 0
       if (mod(iter-1,nprint) == 0) then
         call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
       end if
       if (iflag < 0) go to 300
     end if
!
!  compute the qr factorization of the jacobian.
!
     call qrfac ( m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
!
!  on the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if (iter == 1) then

       if ( mode /= 2 ) then
         diag(1:n) = wa2(1:n)
         do j = 1, n
           if (wa2(j) == 0.0) diag(j) = 1.0
         end do
       end if
!
!  on the first iteration, calculate the norm of the scaled x
!  and initialize the step bound delta.
!
       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm(n,wa3)
       delta = factor*xnorm
       if (delta == 0.0) delta = factor
     end if
!
!  form (q transpose)*fvec and store the first n components in qtf.
!
     wa4(1:m) = fvec(1:m)

     do j = 1, n

       if ( fjac(j,j) /= 0.0 ) then
         sum = 0.0
         do i = j, m
           sum = sum + fjac(i,j)*wa4(i)
         end do
         temp = -sum/fjac(j,j)
         wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
       end if

       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)

     end do
!
!  compute the norm of the scaled gradient.
!
     gnorm = 0.0

     if ( fnorm /= 0.0 ) then

       do j = 1, n
         l = ipvt(j)
         if (wa2(l) /= 0.0) then
           sum = 0.0
           do i = 1, j
             sum = sum + fjac(i,j)*(qtf(i)/fnorm)
           end do
           gnorm = max ( gnorm,abs(sum/wa2(l)))
         end if
       end do

     end if
!
!  test for convergence of the gradient norm.
!
     if ( gnorm <= gtol ) then
       info = 4
       go to 300
     end if
!
!  rescale if necessary.
!
     if (mode /= 2) then
       do j = 1, n
         diag(j) = max ( diag(j),wa2(j))
       end do
     end if
!
!  beginning of the inner loop.
!
  200    continue
!
!  determine the Levenberg-Marquardt parameter.
!
        call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  store the direction p and x + p. calculate the norm of p.
!
        wa1(1:n) = - wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm(n,wa3)
!
!  on the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  evaluate the function at x + p and calculate its norm.
!
        iflag = 1
        call fcn(m,n,wa2,wa4,fjac,ldfjac,iflag)
        nfev = nfev + 1
        if (iflag < 0) go to 300
        fnorm1 = enorm(m,wa4)
!
!  compute the scaled actual reduction.
!
        actred = -1.0
        if (0.1 *fnorm1 < fnorm) actred = 1.0 - (fnorm1/fnorm)**2
!
!  compute the scaled predicted reduction and
!  the scaled directional derivative.
!
        do j = 1, n
          wa3(j) = 0.0
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        end do

        temp1 = enorm(n,wa3)/fnorm
        temp2 = (sqrt(par)*pnorm)/fnorm
        prered = temp1**2 + temp2**2/0.5
        dirder = -(temp1**2 + temp2**2)
!
!  compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0
        if (prered /= 0.0) ratio = actred/prered
!
!  update the step bound.
!
        if (ratio <= 0.25) then
           if (actred >= 0.0) temp = 0.5
           if (actred < 0.0) then
             temp = 0.5*dirder/(dirder + 0.5*actred)
           end if

           if (0.1 *fnorm1 >= fnorm .or. temp < 0.1 ) temp = 0.1
           delta = temp*min ( delta,pnorm/0.1)
           par = par/temp
        else
           if ( par == 0.0 .or. ratio >= 0.75 ) then
             delta = 2.0 * pnorm
             par = 0.5 * par
           end if
        end if
!
!  successful iteration. update x, fvec, and their norms.
!
        if (ratio >= 0.0001 ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm(n,wa2)
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  tests for convergence.
!
        if ( abs(actred) <= ftol .and. prered <= ftol .and. 0.5*ratio <= 1.0) then
          info = 1
        end if

        if (delta <= xtol*xnorm) info = 2
        if (abs(actred) <= ftol .and. prered <= ftol &
          .and. 0.5*ratio <= 1.0 .and. info == 2) then
          info = 3
        end if

        if (info /= 0) go to 300
!
!  tests for termination and stringent tolerances.
!
        if (nfev >= maxfev) info = 5
        if (abs(actred) <= epsmch .and. prered <= epsmch &
          .and. 0.5*ratio <= 1.0) info = 6
        if (delta <= epsmch*xnorm) info = 7
        if (gnorm <= epsmch) info = 8
        if (info /= 0) go to 300
!
!  end of the inner loop. repeat if iteration unsuccessful.
!
        if (ratio < 0.0001 ) go to 200
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

  if ( nprint > 0 ) then
    call fcn(m,n,x,fvec,fjac,ldfjac,iflag)
  end if

  return
end
subroutine lmder1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )
!
!*******************************************************************************
!
!! LMDER1 minimizes M functions in N variables using the Levenberg-Marquardt method.
!
!
!     the purpose of LMDER1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     Levenberg-Marquardt algorithm. this is done by using the more
!     general least-squares solver LMDER. the user must provide a
!     subroutine which calculates the functions and the jacobian.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the jacobian. fcn must
!         be declared in an external statement in the user
!         calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!         integer m,n,ldfjac,iflag
!         real x(n),fvec(m),fjac(ldfjac,n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec. do not alter fjac.
!         if iflag = 2 calculate the jacobian at x and
!         return this matrix in fjac. do not alter fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of LMDER1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached 100*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
  integer ldfjac
  integer m
  integer n
!
  real diag(n)
  real factor
  real fjac(ldfjac,n)
  real ftol
  real fvec(m)
  real gtol
  integer info
  integer ipvt(n)
  integer maxfev
  integer mode
  integer nfev
  integer njev
  integer nprint
  real qtf(n)
  real tol
  real x(n)
  real xtol
!
  external fcn
!
  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( ldfjac < m ) then
    return
  else if ( tol < 0.0 ) then
    return
  end if

  factor = 100.0
  maxfev = 100 * ( n + 1 )
  ftol = tol
  xtol = tol
  gtol = 0.0
  mode = 1
  nprint = 0

  call lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

  if (info == 8) then
    info = 4
  end if

  return
end
subroutine lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
  diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )
!
!*******************************************************************************
!
!! LMDIF minimizes M functions in N variables using the Levenberg-Marquardt method.
!
!
!     the purpose of LMDIF is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the Levenberg-Marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         real x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of LMDIF.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
  integer ldfjac
  integer m
  integer n
!
  real actred
  real delta
  real diag(n)
  real dirder
  real enorm
  real epsfcn
  real epsmch
  real factor
  real fjac(ldfjac,n)
  real fnorm
  real fnorm1
  real ftol
  real fvec(m)
  real gnorm
  real gtol
  integer i
  integer iflag
  integer iter
  integer info
  integer ipvt(n)
  integer j
  integer l
  integer maxfev
  integer mode
  integer nfev
  integer nprint
  real par
  real pnorm
  real prered
  real qtf(n)
  real ratio
  real sum
  real temp
  real temp1
  real temp2
  real wa1(n)
  real wa2(n)
  real wa3(n)
  real wa4(m)
  real x(n)
  real xnorm
  real xtol
!
  external fcn
!
  epsmch = epsilon ( 1.0 )

  info = 0
  iflag = 0
  nfev = 0

  if ( n <= 0 ) then
    go to 300
  else if ( m < n ) then
    go to 300
  else if ( ldfjac < m ) then
    go to 300
  else if ( ftol < 0.0 ) then
    go to 300
  else if ( xtol < 0.0 ) then
    go to 300
  else if ( gtol < 0.0 ) then
    go to 300
  else if ( maxfev <= 0 ) then
    go to 300
  else if ( factor <= 0.0 ) then
    go to 300
  end if

  if (mode == 2) then
    do j = 1, n
      if (diag(j) <= 0.0) go to 300
    end do
  end if
!
!  evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn(m,n,x,fvec,iflag)
  nfev = 1
  if (iflag < 0) go to 300
  fnorm = enorm(m,fvec)
!
!  initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0
  iter = 1
!
!  beginning of the outer loop.
!
   30 continue
!
!  calculate the jacobian matrix.
!
     iflag = 2
     call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
     nfev = nfev + n
     if (iflag < 0) go to 300
!
!  if requested, call fcn to enable printing of iterates.
!
     if ( nprint > 0 ) then
       iflag = 0
       if (mod(iter-1,nprint) == 0) then
         call fcn(m,n,x,fvec,iflag)
       end if
       if (iflag < 0) go to 300
     end if
!
!  compute the qr factorization of the jacobian.
!
     call qrfac ( m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
!
!  on the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if (iter == 1) then

       if (mode /= 2) then
         diag(1:n) = wa2(1:n)
         do j = 1, n
           if (wa2(j) == 0.0) diag(j) = 1.0
         end do
       end if
!
!  on the first iteration, calculate the norm of the scaled x
!  and initialize the step bound delta.
!
       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm(n,wa3)
       delta = factor*xnorm
       if (delta == 0.0) delta = factor
     end if
!
!  form (q transpose)*fvec and store the first n components in qtf.
!
     wa4(1:m) = fvec(1:m)

     do j = 1, n

       if ( fjac(j,j) /= 0.0 ) then
         sum = 0.0
         do i = j, m
           sum = sum + fjac(i,j)*wa4(i)
         end do

         temp = -sum/fjac(j,j)
         wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
       end if

       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)

     end do
!
!  compute the norm of the scaled gradient.
!
     gnorm = 0.0

     if ( fnorm /= 0.0 ) then

       do j = 1, n

         l = ipvt(j)

         if (wa2(l) /= 0.0 ) then
           sum = 0.0
           do i = 1, j
             sum = sum + fjac(i,j)*(qtf(i)/fnorm)
           end do
           gnorm = max ( gnorm,abs ( sum / wa2(l) ) )
         end if

       end do

     end if
!
!  test for convergence of the gradient norm.
!
     if ( gnorm <= gtol ) then
       info = 4
       go to 300
     end if
!
!  rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  beginning of the inner loop.
!
  200    continue
!
!  determine the Levenberg-Marquardt parameter.
!
        call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  store the direction p and x + p. calculate the norm of p.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm(n,wa3)
!
!  on the first iteration, adjust the initial step bound.
!
        if (iter == 1) then
          delta = min ( delta, pnorm )
        end if
!
!  evaluate the function at x + p and calculate its norm.
!
        iflag = 1
        call fcn(m,n,wa2,wa4,iflag)
        nfev = nfev + 1
        if (iflag < 0) go to 300
        fnorm1 = enorm ( m, wa4 )
!
!  compute the scaled actual reduction.
!
        actred = -1.0
        if (0.1 *fnorm1 < fnorm) actred = 1.0 - (fnorm1/fnorm)**2
!
!  compute the scaled predicted reduction and the scaled directional derivative.
!
        do j = 1, n
           wa3(j) = 0.0
           l = ipvt(j)
           temp = wa1(l)
           wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        end do

        temp1 = enorm ( n, wa3 ) / fnorm
        temp2 = ( sqrt ( par ) * pnorm ) / fnorm
        prered = temp1**2 + temp2**2/0.5
        dirder = - ( temp1**2 + temp2**2 )
!
!  compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0
        if (prered /= 0.0) ratio = actred/prered
!
!  update the step bound.
!
        if (ratio <= 0.25 ) then
           if (actred >= 0.0) temp = 0.5
           if (actred < 0.0) then
             temp = 0.5*dirder/(dirder + 0.5*actred)
           end if
           if (0.1 *fnorm1 >= fnorm .or. temp < 0.1 ) temp = 0.1
           delta = temp*min ( delta,pnorm/0.1)
           par = par/temp
        else
           if (par == 0.0 .or. ratio >= 0.75 ) then
             delta = 2.0 * pnorm
             par = 0.5*par
           end if
        end if
!
!  test for successful iteration.
!

!
!  Successful iteration. update x, fvec, and their norms.
!
        if (ratio >= 0.0001 ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm(n,wa2)
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  tests for convergence.
!
        if (abs(actred) <= ftol .and. prered <= ftol &
          .and. 0.5*ratio <= 1.0) then
          info = 1
        end if

        if (delta <= xtol*xnorm) info = 2
        if (abs(actred) <= ftol .and. prered <= ftol &
          .and. 0.5*ratio <= 1.0 .and. info == 2) info = 3
        if (info /= 0) go to 300
!
!  tests for termination and stringent tolerances.
!
        if (nfev >= maxfev) info = 5
        if (abs(actred) <= epsmch .and. prered <= epsmch &
          .and. 0.5*ratio <= 1.0) info = 6
        if (delta <= epsmch*xnorm) info = 7
        if (gnorm <= epsmch) info = 8
        if (info /= 0) go to 300
!
!  end of the inner loop. repeat if iteration unsuccessful.
!
        if (ratio < 0.0001 ) go to 200
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

  if ( nprint > 0 ) then
    call fcn(m,n,x,fvec,iflag)
  end if

  return
end
subroutine lmdif1 ( fcn, m, n, x, fvec, tol, info )
!
!*******************************************************************************
!
!! LMDIF1 minimizes M functions in N variables using the Levenberg-Marquardt method.
!
!
!     the purpose of LMDIF1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     Levenberg-Marquardt algorithm. this is done by using the more
!     general least-squares solver LMDIF. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         real x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of LMDIF1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded 200*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
  integer m
  integer n
!
  real diag(n)
  real epsfcn
  real factor
  real fjac(m,n)
  real ftol
  real fvec(m)
  real gtol
  integer info
  integer ipvt(n)
  integer ldfjac
  integer maxfev
  integer mode
  integer nfev
  integer nprint
  real qtf(n)
  real tol
  real x(n)
  real xtol
!
  external fcn
!
  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( tol < 0.0 ) then
    return
  end if

  factor = 100.0
  maxfev = 200 * ( n + 1 )
  ftol = tol
  xtol = tol
  gtol = 0.0
  epsfcn = 0.0
  mode = 1
  nprint = 0
  ldfjac = m

  call lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
    diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end
subroutine lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )
!
!*******************************************************************************
!
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta,
!     the problem is to determine a value for the parameter
!     par such that if x solves the system
!
!           a*x = b ,     sqrt(par)*d*x = 0 ,
!
!     in the least squares sense, and dxnorm is the euclidean
!     norm of d*x, then either par is zero and
!
!           (dxnorm-delta) <= 0.1*delta ,
!
!     or par is positive and
!
!           abs(dxnorm-delta) <= 0.1*delta .
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then LMPAR expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. on output
!     LMPAR also provides an upper triangular matrix s such that
!
!            t   t                   t
!           p *(a *a + par*d*d)*p = s *s .
!
!     s is employed within LMPAR and may be of separate interest.
!
!     only a few iterations are generally needed for convergence
!     of the algorithm. if, however, the limit of 10 iterations
!     is reached, then the output par will contain the best
!     value obtained so far.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
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
!       par is a nonnegative variable. on input par contains an
!         initial estimate of the Levenberg-Marquardt parameter.
!         on output par contains the final estimate.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
!         for the output par.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
  integer ldr
  integer n
!
  real delta
  real diag(n)
  real dwarf
  real dxnorm
  real enorm
  real gnorm
  real fp
  integer i
  integer ipvt(n)
  integer iter
  integer j
  integer k
  integer l
  integer nsing
  real par
  real parc
  real parl
  real paru
  real qnorm
  real qtb(n)
  real r(ldr,n)
  real sdiag(n)
  real sum
  real temp
  real wa1(n)
  real wa2(n)
  real x(n)
!
!  dwarf is the smallest positive magnitude.
!
  dwarf = tiny ( 1.0 )
!
!  compute and store in x the gauss-newton direction. if the
!  jacobian is rank-deficient, obtain a least squares solution.
!
  nsing = n

  do j = 1, n
     wa1(j) = qtb(j)
     if (r(j,j) == 0.0 .and. nsing == n) nsing = j - 1
     if (nsing < n) wa1(j) = 0.0
  end do

  do k = 1, nsing
     j = nsing - k + 1
     wa1(j) = wa1(j)/r(j,j)
     temp = wa1(j)
     wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  end do

  do j = 1, n
     l = ipvt(j)
     x(l) = wa1(j)
  end do
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the gauss-newton direction.
!
  iter = 0
  wa2(1:n) = diag(1:n) * x(1:n)
  dxnorm = enorm(n,wa2)
  fp = dxnorm - delta

  if (fp <= 0.1 *delta) go to 220
!
!  if the jacobian is not rank deficient, the newton
!  step provides a lower bound, parl, for the zero of
!  the function. otherwise set this bound to zero.
!
  parl = 0.0

  if ( nsing >= n ) then

    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l)*(wa2(l)/dxnorm)
    end do

    do j = 1, n
      sum = 0.0
      do i = 1, j-1
        sum = sum + r(i,j)*wa1(i)
      end do
      wa1(j) = (wa1(j) - sum)/r(j,j)
    end do

    temp = enorm ( n, wa1 )
    parl = (( fp / delta ) / temp ) / temp

  end if
!
!  calculate an upper bound, paru, for the zero of the function.
!
  do j = 1, n
     sum = 0.0
     do i = 1, j
        sum = sum + r(i,j)*qtb(i)
     end do
     l = ipvt(j)
     wa1(j) = sum/diag(l)
  end do

  gnorm = enorm(n,wa1)
  paru = gnorm/delta
  if (paru == 0.0) paru = dwarf / min ( delta,0.1 )
!
!  if the input par lies outside of the interval (parl,paru),
!  set par to the closer endpoint.
!
  par = max ( par,parl)
  par = min ( par,paru)
  if (par == 0.0) par = gnorm/dxnorm
!
!  beginning of an iteration.
!
  150 continue

     iter = iter + 1
!
!  evaluate the function at the current value of par.
!
     if (par == 0.0) par = max ( dwarf, 0.001 * paru )

     wa1(1:n) = sqrt ( par ) * diag(1:n)

     call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

     wa2(1:n) = diag(1:n)*x(1:n)
     dxnorm = enorm(n,wa2)
     temp = fp
     fp = dxnorm - delta
!
!  if the function is small enough, accept the current value
!  of par. also test for the exceptional cases where parl
!  is zero or the number of iterations has reached 10.
!
     if (abs(fp) <= 0.1 *delta .or. parl == 0.0 .and. fp <= temp &
       .and. temp < 0.0 .or. iter == 10) go to 220
!
!  compute the newton correction.
!
     do j = 1, n
        l = ipvt(j)
        wa1(j) = diag(l)*(wa2(l)/dxnorm)
     end do

     do j = 1, n
        wa1(j) = wa1(j) / sdiag(j)
        temp = wa1(j)
        wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
     end do

     temp = enorm(n,wa1)
     parc = ((fp/delta)/temp)/temp
!
!  depending on the sign of the function, update parl or paru.
!
     if ( fp > 0.0 ) then
       parl = max ( parl,par)
     else if (fp < 0.0) then
       paru = min ( paru,par)
     end if
!
!  compute an improved estimate for par.
!
     par = max ( parl,par+parc)
!
!  end of an iteration.
!
     go to 150
  220 continue
!
!  termination.
!
  if ( iter == 0 ) then
    par = 0.0
  end if

  return
end
subroutine lmstr ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
  diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )
!
!*******************************************************************************
!
!! LMSTR minimizes M functions in N variables using the Levenberg-Marquardt method.
!
!
!     the purpose of LMSTR is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the Levenberg-Marquardt algorithm which uses minimal storage.
!     the user must provide a subroutine which calculates the
!     functions and the rows of the jacobian.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the rows of the jacobian.
!         fcn must be declared in an external statement in the
!         user calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjrow,iflag)
!         integer m,n,iflag
!         real x(n),fvec(m),fjrow(n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec.
!         if iflag = i calculate the (i-1)-st row of the
!         jacobian at x and return this vector in fjrow.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of LMSTR.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array. the upper triangle of fjac
!         contains an upper triangular matrix r such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower triangular
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn with iflag = 1
!         has reached maxfev.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn with iflag = 1.
!
!       njev is an integer output variable set to the number of
!         calls to fcn with iflag = 2.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
  integer ldfjac
  integer m
  integer n
!
  real actred
  real delta
  real diag(n)
  real dirder
  real enorm
  real epsmch
  real factor
  real fjac(ldfjac,n)
  real fnorm
  real fnorm1
  real ftol
  real fvec(m)
  real gnorm
  real gtol
  integer i
  integer iflag
  integer info
  integer ipvt(n)
  integer iter
  integer j
  integer l
  integer maxfev
  integer mode
  integer nfev
  integer njev
  integer nprint
  real par
  real pnorm
  real prered
  real qtf(n)
  real ratio
  logical sing
  real sum
  real temp
  real temp1
  real temp2
  real wa1(n)
  real wa2(n)
  real wa3(n)
  real wa4(m)
  real x(n)
  real xnorm
  real xtol
!
  external fcn
!
  epsmch = epsilon ( 1.0 )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  check the input parameters for errors.
!
  if ( n <= 0 ) then
    go to 340
  else if ( m < n ) then
    go to 340 
  else if ( ldfjac < n ) then
    go to 340
  else if ( ftol < 0.0 ) then
    go to 340
  else if ( xtol < 0.0 ) then
    go to 340
  else if ( gtol < 0.0 ) then
    go to 340
  else if ( maxfev <= 0 ) then
    go to 340
  else if ( factor <= 0.0 ) then
    go to 340
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if (diag(j) <= 0.0) go to 340
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  call fcn(m,n,x,fvec,wa3,iflag)
  nfev = 1
  if (iflag < 0) go to 340
  fnorm = enorm(m,fvec)
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0
  iter = 1
!
!  beginning of the outer loop.
!
   30 continue
!
!  If requested, call fcn to enable printing of iterates.
!
     if (nprint > 0) then
       iflag = 0
       if (mod(iter-1,nprint) == 0) then
         call fcn(m,n,x,fvec,wa3,iflag)
       end if
       if (iflag < 0) go to 340
     end if
!
!  Compute the qr factorization of the jacobian matrix calculated one row 
!  at a time, while simultaneously forming (q transpose)*fvec and storing 
!  the first n components in qtf.
!
     qtf(1:n) = 0.0
     fjac(1:n,1:n) = 0.0
     iflag = 2

     do i = 1, m
        call fcn(m,n,x,fvec,wa3,iflag)
        if (iflag < 0) go to 340
        temp = fvec(i)
        call rwupdt(n,fjac,ldfjac,wa3,qtf,temp,wa1,wa2)
        iflag = iflag + 1
     end do

     njev = njev + 1
!
!  If the jacobian is rank deficient, call QRFAC to
!  reorder its columns and update the components of qtf.
!
     sing = .false.
     do j = 1, n
        if (fjac(j,j) == 0.0) sing = .true.
        ipvt(j) = j
        wa2(j) = enorm(j,fjac(1,j))
     end do

     if ( sing ) then

       call qrfac ( n, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )

       do j = 1, n

         if ( fjac(j,j) /= 0.0 ) then

           sum = 0.0
           do i = j, n
             sum = sum + fjac(i,j)*qtf(i)
           end do

           temp = - sum / fjac(j,j)
           qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp

         end if

         fjac(j,j) = wa1(j)

       end do

     end if
!
!  on the first iteration 
!    if mode is 1, 
!      scale according to the norms of the columns of the initial jacobian.
!    calculate the norm of the scaled x,
!    initialize the step bound delta.
!
     if ( iter == 1 ) then

       if (mode /= 2) then

         diag(1:n) = wa2(1:n)
         do j = 1, n
           if (wa2(j) == 0.0) diag(j) = 1.0
         end do

       end if

       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm(n,wa3)
       delta = factor*xnorm
       if (delta == 0.0) delta = factor

     end if
!
!  compute the norm of the scaled gradient.
!
     gnorm = 0.0

     if ( fnorm /= 0.0 ) then

       do j = 1, n
         l = ipvt(j)
         if ( wa2(l) /= 0.0 ) then
           sum = 0.0
           do i = 1, j
             sum = sum + fjac(i,j) * ( qtf(i) / fnorm )
           end do
           gnorm = max ( gnorm, abs ( sum / wa2(l) ) )
         end if
       end do

     end if
!
!  test for convergence of the gradient norm.
!
     if (gnorm <= gtol) then
       info = 4
       go to 340
     end if
!
!  rescale if necessary.
!
     if ( mode /= 2) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  beginning of the inner loop.
!
  240    continue
!
!  determine the Levenberg-Marquardt parameter.
!
        call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  store the direction p and x + p. calculate the norm of p.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)
        pnorm = enorm(n,wa3)
!
!  on the first iteration, adjust the initial step bound.
!
        if (iter == 1) delta = min ( delta,pnorm)
!
!  evaluate the function at x + p and calculate its norm.
!
        iflag = 1
        call fcn(m,n,wa2,wa4,wa3,iflag)
        nfev = nfev + 1
        if (iflag < 0) go to 340
        fnorm1 = enorm(m,wa4)
!
!  compute the scaled actual reduction.
!
        actred = -1.0
        if (0.1 *fnorm1 < fnorm) actred = 1.0 - (fnorm1/fnorm)**2
!
!  compute the scaled predicted reduction and
!  the scaled directional derivative.
!
        do j = 1, n
           wa3(j) = 0.0
           l = ipvt(j)
           temp = wa1(l)
           wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        end do

        temp1 = enorm(n,wa3)/fnorm
        temp2 = (sqrt(par)*pnorm)/fnorm
        prered = temp1**2 + temp2**2/0.5
        dirder = -(temp1**2 + temp2**2)
!
!  compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0
        if (prered /= 0.0) ratio = actred/prered
!
!  update the step bound.
!
        if (ratio <= 0.25 ) then
           if (actred >= 0.0) temp = 0.5
           if (actred < 0.0) then
             temp = 0.5*dirder/(dirder + 0.5*actred)
           end if
           if (0.1 *fnorm1 >= fnorm .or. temp < 0.1 ) temp = 0.1
           delta = temp*min ( delta,pnorm/0.1 )
           par = par/temp
        else
           if ( par == 0.0 .or. ratio >= 0.75 ) then
             delta = pnorm/0.5
             par = 0.5*par
           end if
        end if
!
!  test for successful iteration.
!
        if (ratio >= 0.0001 ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm(n,wa2)
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  tests for convergence, termination and stringent tolerances.
!
        if (abs(actred) <= ftol .and. prered <= ftol &
          .and. 0.5*ratio <= 1.0) then
          info = 1
        end if

        if ( delta <= xtol*xnorm ) then
          info = 2
        end if

        if (abs(actred) <= ftol .and. prered <= ftol &
          .and. 0.5*ratio <= 1.0 .and. info == 2) then
          info = 3
        end if

        if ( info /= 0 ) then
          go to 340
        end if

        if (nfev >= maxfev) info = 5
        if (abs(actred) <= epsmch .and. prered <= epsmch &
          .and. 0.5*ratio <= 1.0) info = 6
        if (delta <= epsmch*xnorm) info = 7
        if (gnorm <= epsmch) info = 8
        if (info /= 0) go to 340
!
!  end of the inner loop. repeat if iteration unsuccessful.
!
        if (ratio < 0.0001 ) go to 240
!
!  end of the outer loop.
!
     go to 30

  340 continue
!
!  termination, either normal or user imposed.
!
  if (iflag < 0) info = iflag
  iflag = 0

  if ( nprint > 0 ) then
    call fcn(m,n,x,fvec,wa3,iflag)
  end if

  return
end
subroutine lmstr1 ( fcn, m, n, x, fvec, fjac, ldfjac, tol, info )
!
!*******************************************************************************
!
!! LMSTR1 minimizes M functions in N variables using the Levenberg-Marquardt method.
!
!
!     the purpose of LMSTR1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the Levenberg-Marquardt algorithm which uses minimal storage.
!     this is done by using the more general least-squares solver
!     LMSTR. the user must provide a subroutine which calculates
!     the functions and the rows of the jacobian.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions and the rows of the jacobian.
!         fcn must be declared in an external statement in the
!         user calling program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,fjrow,iflag)
!         integer m,n,iflag
!         real x(n),fvec(m),fjrow(n)
!         ----------
!         if iflag = 1 calculate the functions at x and
!         return this vector in fvec.
!         if iflag = i calculate the (i-1)-st row of the
!         jacobian at x and return this vector in fjrow.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of LMSTR1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       fjac is an output n by n array. the upper triangle of fjac
!         contains an upper triangular matrix r such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower triangular
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn with iflag = 1 has
!                   reached 100*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
  integer ldfjac
  integer m
  integer n
!
  real diag(n)
  real factor
  real fjac(ldfjac,n)
  real ftol
  real fvec(m)
  real gtol
  integer info
  integer ipvt(n)
  integer maxfev
  integer mode
  integer nfev
  integer njev
  integer nprint
  real qtf(n)
  real tol
  real x(n)
  real xtol
!
  external fcn
!
  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( ldfjac < n ) then
    return
  else if ( tol < 0.0 ) then
    return
  end if

  factor = 100.0
  maxfev = 100*(n + 1)
  ftol = tol
  xtol = tol
  gtol = 0.0
  mode = 1
  nprint = 0

  call lmstr ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
    diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf )

  if ( info == 8 ) then
    info = 4
  end if

  return
end
subroutine qform ( m, n, q, ldq )
!
!*******************************************************************************
!
!! QFORM constructs an orthogonal matrix Q from its factored form.
!
!
!     this subroutine proceeds from the computed qr factorization of
!     an m by n matrix a to accumulate the m by m orthogonal matrix
!     q from its factored form.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
  integer ldq
  integer m
!
  integer i
  integer j
  integer k
  integer l
  integer minmn
  integer n
  real q(ldq,m)
  real sum
  real temp
  real wa(m)
!
!  Zero out upper triangle of q in the first min(m,n) columns.
!
  minmn = min ( m, n )

  do j = 2, minmn
    q(1:j-1,j) = 0.0
  end do
!
!  Initialize remaining columns to those of the identity matrix.
!
  do j = n+1, m
    q(1:m,j) = 0.0
    q(j,j) = 1.0
  end do
!
!  Accumulate Q from its factored form.
!
  do l = 1, minmn

     k = minmn - l + 1

     wa(k:m) = q(k:m,k)
     q(k:m,k) = 0.0
     q(k,k) = 1.0

     if ( wa(k) /= 0.0 ) then

      do j = k, m
        sum = 0.0
        do i = k, m
          sum = sum + q(i,j)*wa(i)
        end do
        temp = sum / wa(k)
        q(k:m,j) = q(k:m,j) - temp * wa(k:m)
      end do

    end if

  end do

  return
end
subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm ) 
!
!*******************************************************************************
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!
!     this subroutine uses Householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the Householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           t
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
!       pivot is a logical input variable. if pivot is set true,
!         then column pivoting is enforced. if pivot is set false,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.
!
!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.
!
!       rdiag is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.
!
  integer lda
  integer lipvt
  integer n
!
  real a(lda,n)
  real acnorm(n)
  real ajnorm
  real enorm
  real epsmch
  integer i
  integer ipvt(lipvt)
  integer j
  integer k
  integer kmax
  integer m
  integer minmn
  logical pivot
  real rdiag(n)
  real sum
  real temp
  real wa(n)
!
  epsmch = epsilon ( 1.0 )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, a(1,j) )
  end do

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  if ( pivot ) then
    do j = 1, n
      ipvt(j) = j
    end do
  end if
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pivot ) then

      kmax = j

      do k = j, n
        if (rdiag(k) > rdiag(kmax)) kmax = k
      end do

      if ( kmax /= j ) then

        do i = 1, m
          call r_swap ( a(i,j), a(i,kmax) )
        end do

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)
        k = ipvt(j)
        ipvt(j) = ipvt(kmax)
        ipvt(kmax) = k

      end if

    end if
!
!  Compute the Householder transformation to reduce the
!  j-th column of a to a multiple of the j-th unit vector.
!
     ajnorm = enorm(m-j+1,a(j,j))

     if (ajnorm /= 0.0) then

     if (a(j,j) < 0.0) ajnorm = -ajnorm
     a(j:m,j) = a(j:m,j) / ajnorm
     a(j,j) = a(j,j) + 1.0
!
!  Apply the transformation to the remaining columns and update the norms.
!
     do k = j+1, n

        sum = 0.0
        do i = j, m
          sum = sum + a(i,j)*a(i,k)
        end do

        temp = sum / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        if ( pivot .and. rdiag(k) /= 0.0 ) then

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k)*sqrt(max ( 0.0,1.0-temp**2))

          if (0.05*(rdiag(k)/wa(k))**2 <= epsmch) then
            rdiag(k) = enorm(m-j,a(j+1,k))
            wa(k) = rdiag(k)
          end if

        end if

     end do

     end if

     rdiag(j) = -ajnorm

    end do

  return
end
subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )
!
!*******************************************************************************
!
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!
!     given an m by n matrix a, an n by n diagonal matrix d,
!     and an m-vector b, the problem is to determine an x which
!     solves the system
!
!           a*x = b ,     d*x = 0 ,
!
!     in the least squares sense.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then QRSOLV expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. the system
!     a*x = b, d*x = 0, is then equivalent to
!
!                  t       t
!           r*z = q *b ,  p *d*p*z = 0 ,
!
!     where x = p*z. if this system does not have full rank,
!     then a least squares solution is obtained. on output QRSOLV
!     also provides an upper triangular matrix s such that
!
!            t   t               t
!           p *(a *a + d*d)*p = s *s .
!
!     s is computed within QRSOLV and may be of separate interest.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, d*x = 0.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
  integer ldr
  integer n
!
  real c
  real cotan
  real diag(n)
  integer i
  integer ipvt(n)
  integer j
  integer k
  integer l
  integer nsing
  real qtb(n)
  real qtbpj
  real r(ldr,n)
  real s
  real sdiag(n)
  real sum
  real t
  real temp
  real wa(n)
  real x(n)
!
!  copy r and (q transpose)*b to preserve input and initialize s.
!  in particular, save the diagonal elements of r in x.
!
  do j = 1, n
     r(j:n,j) = r(j,j:n)
     x(j) = r(j,j)
  end do

  wa(1:n) = qtb(1:n)
!
!  eliminate the diagonal matrix d using a Givens rotation.
!
  do j = 1, n
!
!  Prepare the row of d to be eliminated, locating the
!  diagonal element using p from the qr factorization.
!
    l = ipvt(j)

    if ( diag(l) /= 0.0 ) then

      sdiag(j:n) = 0.0
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of d
!  modify only a single element of (q transpose)*b
!  beyond the first n, which is initially zero.
!
      qtbpj = 0.0

      do k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of d.
!
        if ( sdiag(k) /= 0.0 ) then

          if ( abs ( r(k,k) ) < abs ( sdiag(k) ) ) then
            cotan = r(k,k) / sdiag(k)
            s = 0.5 / sqrt ( 0.25 + 0.25 * cotan**2 )
            c = s * cotan
          else
            t = sdiag(k) / r(k,k)
            c = 0.5 / sqrt ( 0.25 + 0.25 * t**2 )
            s = c * t
          end if
!
!  Compute the modified diagonal element of r and
!  the modified element of ((q transpose)*b,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          do i = k+1, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          end do

        end if

      end do

    end if
!
!  Store the diagonal element of s and restore
!  the corresponding diagonal element of r.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  end do
!
!  Solve the triangular system for z. if the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  do j = 1, n
    if (sdiag(j) == 0.0 .and. nsing == n) nsing = j - 1
    if (nsing < n) wa(j) = 0.0
  end do

  do k = 1, nsing
    j = nsing - k + 1
    sum = 0.0
    do i = j+1, nsing
      sum = sum + r(i,j)*wa(i)
    end do
    wa(j) = (wa(j) - sum)/sdiag(j)
  end do
!
!  Permute the components of z back to components of x.
!
  do j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
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
subroutine rwupdt ( n, r, ldr, w, b, alpha, c, s )
!
!*******************************************************************************
!
!! RWUPDT computes the decomposition of a triangular matrix augmented by one row.
!
!
!     given an n by n upper triangular matrix r, this subroutine
!     computes the qr decomposition of the matrix formed when a row
!     is added to r. if the row is specified by the vector w, then
!     rwupdt determines an orthogonal matrix q such that when the
!     n+1 by n matrix composed of r augmented by w is premultiplied
!     by (q transpose), the resulting matrix is upper trapezoidal.
!     the matrix (q transpose) is the product of n transformations
!
!           g(n)*g(n-1)* ... *g(1)
!
!     where g(i) is a Givens rotation in the (i,n+1) plane which
!     eliminates elements in the (n+1)-st plane. rwupdt also
!     computes the product (q transpose)*c where c is the
!     (n+1)-vector (b,alpha). q itself is not accumulated, rather
!     the information to recover the g rotations is supplied.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
!
!  Parameters:
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the upper triangular part of
!         r must contain the matrix to be updated. on output r
!         contains the updated triangular matrix.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       w is an input array of length n which must contain the row
!         vector to be added to r.
!
!       b is an array of length n. on input b must contain the
!         first n elements of the vector c. on output b contains
!         the first n elements of the vector (q transpose)*c.
!
!       alpha is a variable. on input alpha must contain the
!         (n+1)-st element of the vector c. on output alpha contains
!         the (n+1)-st element of the vector (q transpose)*c.
!
!       c is an output array of length n which contains the
!         cosines of the transforming Givens rotations.
!
!       s is an output array of length n which contains the
!         sines of the transforming Givens rotations.
!
  integer ldr
  integer n
!
  real alpha
  real b(n)
  real c(n)
  real cotan
  integer i
  integer j
  real r(ldr,n)
  real rowj
  real s(n)
  real tan
  real temp
  real w(n)
!
  do j = 1, n

    rowj = w(j)
!
!  apply the previous transformations to r(i,j), i=1,2,...,j-1, and to w(j).
!
    do i = 1, j-1
      temp =   c(i) * r(i,j) + s(i) * rowj
      rowj = - s(i) * r(i,j) + c(i) * rowj
      r(i,j) = temp
    end do
!
!  determine a Givens rotation which eliminates w(j).
!
    c(j) = 1.0
    s(j) = 0.0

    if ( rowj /= 0.0 ) then

      if ( abs ( r(j,j) ) < abs ( rowj ) ) then
        cotan = r(j,j) / rowj
        s(j) = 0.5 / sqrt(0.25+0.25*cotan**2)
        c(j) = s(j) * cotan
      else
        tan = rowj / r(j,j)
        c(j) = 0.5 / sqrt(0.25+0.25*tan**2)
        s(j) = c(j) * tan
      end if
!
!  apply the current transformation to r(j,j), b(j), and alpha.
!
      r(j,j) =  c(j) * r(j,j) + s(j) * rowj
      temp =    c(j) * b(j)   + s(j) * alpha
      alpha = - s(j) * b(j)   + c(j) * alpha
      b(j) = temp
   
    end if

  end do

  return
end
subroutine r1mpyq ( m, n, a, lda, v, w )
!
!*******************************************************************************
!
!! R1MPYQ computes A*Q, where Q is the product of Householder transformations.
!
!
!     given an m by n matrix a, this subroutine computes a*q where
!     q is the product of 2*(n - 1) transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     and gv(i), gw(i) are Givens rotations in the (i,n) plane which
!     eliminate elements in the i-th and n-th planes, respectively.
!     q itself is not given, rather the information to recover the
!     gv, gw rotations is supplied.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
!         information necessary to recover the Givens rotation gv(i)
!         described above.
!
!       w is an input array of length n. w(i) must contain the
!         information necessary to recover the Givens rotation gw(i)
!         described above.
!
  integer lda
  integer m
  integer n
!
  real a(lda,n)
  real c
  integer i
  integer j
  real s
  real temp
  real v(n)
  real w(n)
!
!  apply the first set of Givens rotations to a.
!
  do j = n-1, 1, -1

     if ( abs ( v(j) ) > 1.0 ) then
       c = 1.0 / v(j)
       s = sqrt ( 1.0 - c**2 )
     else
       s = v(j)
       c = sqrt ( 1.0 - s**2 )
     end if

     do i = 1, m
        temp =   c * a(i,j) - s * a(i,n)
        a(i,n) = s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do
!
!  apply the second set of Givens rotations to a.
!
  do j = 1, n-1

     if ( abs ( w(j) ) > 1.0 ) then
       c = 1.0 / w(j)
       s = sqrt ( 1.0 - c**2 )
     else
       s = w(j)
       c = sqrt ( 1.0 - s**2 )
     end if

     do i = 1, m
        temp =     c * a(i,j) + s * a(i,n)
        a(i,n) = - s * a(i,j) + c * a(i,n)
        a(i,j) = temp
     end do

  end do

  return
end
subroutine r1updt ( m, n, s, ls, u, v, w, sing )
!
!*******************************************************************************
!
!! R1UPDT re-triangularizes a matrix after a rank one update.
!
!
!     given an m by n lower trapezoidal matrix s, an m-vector u,
!     and an n-vector v, the problem is to determine an
!     orthogonal matrix q such that
!
!                   t
!           (s + u*v )*q
!
!     is again lower trapezoidal.
!
!     this subroutine determines q as the product of 2*(n - 1)
!     transformations
!
!           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
!     where gv(i), gw(i) are Givens rotations in the (i,n) plane
!     which eliminate elements in the i-th and n-th planes,
!     respectively. q itself is not accumulated, rather the
!     information to recover the gv, gw rotations is returned.
!
!  Reference:
!
!    More, Garbow and Hillstrom,
!    User Guide for MINPACK-1
!    Argonne National Laboratory,
!    Argonne, Illinois.
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
!         recover the Givens rotation gv(i) described above.
!
!       w is an output array of length m. w(i) contains information
!         necessary to recover the Givens rotation gw(i) described
!         above.
!
!       sing is a logical output variable. sing is set true if any
!         of the diagonal elements of the output s are zero. otherwise
!         sing is set false.
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
!  giant is the largest magnitude.
!
  giant = huge ( 1.0 )
!
!  initialize the diagonal element pointer.
!
  jj = (n*(2*m - n + 1))/2 - (m - n)
!
!  move the nontrivial part of the last column of s into w.
!
  l = jj
  do i = n, m
    w(i) = s(l)
    l = l + 1
  end do
!
!  rotate the vector v into a multiple of the n-th unit vector
!  in such a way that a spike is introduced into w.
!
  do j = n-1, 1, -1

    jj = jj - (m - j + 1)
    w(j) = 0.0

    if (v(j) /= 0.0) then
!
!  determine a Givens rotation which eliminates the
!  j-th element of v.
!
      if ( abs ( v(n) ) < abs ( v(j) ) ) then
        cotan = v(n) / v(j)
        sin = 0.5 / sqrt ( 0.25 + 0.25 * cotan**2 )
        cos = sin * cotan
        tau = 1.0
        if ( abs ( cos ) * giant > 1.0 ) tau = 1.0 / cos
      else
        tan = v(j) / v(n)
        cos = 0.5 / sqrt ( 0.25 + 0.25 * tan**2 )
        sin = cos * tan
        tau = sin
      end if
!
!  apply the transformation to v and store the information
!  necessary to recover the Givens rotation.
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
!  add the spike from the rank 1 update to w.
!
   w(1:m) = w(1:m) + v(n) * u(1:m)
!
!  eliminate the spike.
!
  sing = .false.

  do j = 1, n-1

    if ( w(j) /= 0.0 ) then
!
!  determine a Givens rotation which eliminates the
!  j-th element of the spike.
!
      if ( abs ( s(jj) ) < abs ( w (j) ) ) then
        cotan = s(jj)/w(j)
        sin = 0.5 /sqrt(0.25+0.25*cotan**2)
        cos = sin*cotan
        tau = 1.0
        if (abs(cos)*giant > 1.0) tau = 1.0/cos
      else
        tan = w(j)/s(jj)
        cos = 0.5 /sqrt(0.25+0.25*tan**2)
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
!  store the information necessary to recover the Givens rotation.
!
      w(j) = tau
  
    end if
!
!  test for zero diagonal elements in the output s.
!
    if ( s(jj) == 0.0 ) then
      sing = .true.
    end if

    jj = jj + (m - j + 1)
  
  end do
!
!  Move W back into the last column of the output S.
!
  l = jj
  do i = n, m
    s(l) = w(i)
    l = l + 1
  end do

  if ( s(jj) == 0.0 ) then
    sing = .true.
  end if

  return
end
