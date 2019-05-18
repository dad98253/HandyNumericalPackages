!  nspcg.f90  22 June 2000
!
subroutine adinfn (nn,ndim,maxnzz,jcoef,coef,nstore,ainf,wksp)
!
!*******************************************************************************
!
!! ADINFN computes an upper bound on the spectral radius of inverse(D)*A.
!
!
!  Parameters:
!
!         n       order of system (= nn)
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array (= maxnzz)
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         nstore  matrix storage mode
!                  = 2   symmetric diagonal format
!                  = 3   nonsymmetric diagonal format
!         ainf    upper bound estimate upon output
!         wksp    workspace vector of length n
!
  integer ndim  
!
  real ainf
  real coef(ndim,1)
  integer i
  integer j
  integer jcoef(2)
  integer jd
  real wksp(1)
!
  n = nn
  maxnz = maxnzz

  if ( ainf > 0.0 ) then
    return
  end if

  do i = 1, n
    wksp(i) = coef(i,1)
  end do

  do 25 jd = 1, maxnz

    do j = 1, maxnz

      if ( jcoef(j) == jd ) then

        do i = 1,n
          wksp(i) = wksp(i) - abs (coef(i,j))
        end do

        if (nstore == 3) go to 25

        do i = 1, n - jd
          wksp(i+jd) = wksp(i+jd) - abs (coef(i,j))
        end do

        go to 25

      end if

    end do

    go to 30

 25   continue

 30   continue

  if (nstore == 2) go to 50

  do 45 jd = 1,maxnz
    do j = 1,maxnz
      if ( jcoef(j) == -jd ) then
        do i = 1,n
          wksp(i) = wksp(i) - abs (coef(i,j))
        end do
        go to 45
      end if
    end do
    go to 50
 45   continue
!
!  factor.
!
 50   continue

  t1 = vmin (n,wksp)
  if (t1 <= 0.0) t1 = 1.0

  call ainfn (n,ndim,maxnz,jcoef,coef,nstore,ainf,wksp)

  ainf = ainf / t1

  return
end
subroutine adjust (n,ndim,maxnzz,jcoef,key)
!
!*******************************************************************************
!
!! ADJUST makes adjustments to the JCOEF array.
!
!
!  Parameters:
!
!          n      dimension of the matrix.
!          ndim   row dimension of jcoef array in defining routine
!          maxnz  number of columns in jcoef array
!          jcoef  integer matrix representation array
!          key    indicator flag
!                  = 1   remove zeros from jcoef array
!                  = 2   restore zeros to jcoef array
!
  integer ndim  
!
  integer jcoef(ndim,1)
!
  maxnz = maxnzz
  if (maxnz < 2) return
  if (key == 2) go to 20
!
!  change zero elements of jcoef array.
!
  do j = 2,maxnz
    do i = 1,n
      if (jcoef(i,j) <= 0) jcoef(i,j) = i
    end do
  end do
 
  return
!
!  put original zeros back in jcoef array.
!
 20   do 30 j = 2,maxnz
     do 25 i = 1,n
 25      if (jcoef(i,j) == i) jcoef(i,j) = 0
 30   continue
  return
end
subroutine ainfn (nn,ndim,maxnzz,jcoef,coef,nstore,ainf,wksp)
!
!*******************************************************************************
!
!! AINFN calculates the infinity norm of a matrix.
!
!
!  Parameters:
!
!         n       order of system (= nn)
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array (= maxnzz)
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         nstore  matrix storage mode
!                  = 1   Purdue format
!                  = 2   symmetric diagonal format
!                  = 3   nonsymmetric diagonal format
!                  = 4   symmetric sparse format
!                  = 5   nonsymmetric sparse format
!         ainf    the infinity norm of the matrix, //a//, upon
!                  output
!         wksp    workspace vector of length n
!
  integer ndim
!
  dimension coef(ndim,1)
  integer jcoef(ndim,2)
  dimension wksp(1)
!
  n = nn
  maxnz = maxnzz

  if (ainf > 0.0) return

  go to (10,30,55,75,75), nstore
!
!  ellpack data structure.
!
 10   continue

  do i = 1,n
    wksp(i) = abs (coef(i,1))
  end do

  do j = 2,maxnz
    do i = 1,n
      wksp(i) = wksp(i) + abs (coef(i,j))
    end do
  end do

  go to 995
!
!  symmetric diagonal data structure.
!
 30   continue

  do i = 1,n
    wksp(i) = abs (coef(i,1))
  end do

  do j = 2,maxnz

    ind = jcoef(j,1)
    len = n - ind

    do i = 1,len
      wksp(i) = wksp(i) + abs (coef(i,j))
    end do

    do i = 1,len
      wksp(i+ind) = wksp(i+ind) + abs (coef(i,j))
    end do

  end do

  go to 995
!
!  nonsymmetric diagonal data structure.
!
 55   continue

  do i = 1,n
    wksp(i) = abs (coef(i,1))
  end do

  do j = 2,maxnz
    ind = jcoef(j,1)
    len = n - iabs(ind)
    ist1 = max(1,1 - ind)
    ist2 = min(n,n - ind)
    do i = ist1,ist2
      wksp(i) = wksp(i) + abs (coef(i,j))
    end do
  end do

  go to 995
!
!  sparse structure.
!
 75   continue

  do i = 1,n
    wksp(i) = abs (coef(i,1))
  end do

  do k = n+1,maxnz
    wksp(jcoef(k,1)) = wksp(jcoef(k,1)) + abs (coef(k,1))
  end do

  if (nstore == 5) go to 995

  do k = n+1,maxnz
    wksp(jcoef(k,2)) = wksp(jcoef(k,2)) + abs (coef(k,1))
  end do
!
!  determine ainf = max (wksp(i)).
!
 995  continue

  ainf = vmax (n,wksp)

  return
end
subroutine basic (suba,subat,subql,subqlt,subqr,subqrt,subadp, &
  coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BASIC is the user interface to the basic (unaccelerated) iterative method, with preconditioning.  
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2)
  dimension wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call basicw ( suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs, &
    wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine basicw ( suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs,wk, &
  nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BASICW runs the basic (unaccelerated) iterative method, with preconditioning.  
!
!
!  that is, it applies the fixed point method
!  to the preconditioned system.
!  two-sided preconditioning is efficiently implemented.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2)
  dimension wfac(1), jwfac(1)
  logical iql, iqr
  external suba, subql, subqr
  dimension iparm(30), rparm(30)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, &
                    alphao, gamma, sigma, rr, rho, dkq, dkm1, &
                    ff, rqmin, rqmax, stptst, udnm, ubarnm, &
                    bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, &
                    udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav, &
                    rdot, rzdot, rztdot, zdot, zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! preliminary calculations
!
  iacel = 0
  ier = 0
  nwusd = 0
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 996
  if (level >= 2) write (nout,496)
496   format (' basic')
! use knowledge about spectrum to optimally extrapolate.
  extrap = (emax+emin)/2.0
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
!
! initialize the stopping test.
!
  call inithv (0)
  zthave = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr, coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 735
!
! bust up workspace.
!
  izt = 1
  iv1 = izt + n
  iwfree = iv1 + n
  if (iqlr == 0) iwfree = iv1
  nwusd = max(nwusd,iwfree-1)
!
! check the memory usage.
!
  if (nwusd > nw) go to 999
!
! do preliminary calculations.
!
  in = 0
  is = 0
  go to (151,152,153,154),iqlr + 1
!
 151  call suba (coef,jcoef,wfac,jwfac,n,u,wk(izt))
  call vexopy (n,wk(izt),rhs,wk(izt),2)
  go to 10
!
 152  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(izt))
  go to 10
!
 153  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subqr (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(izt))
  go to 10
!
 154  call suba (coef,jcoef,wfac,jwfac,n,u,wk(izt))
  call vexopy (n,wk(izt),rhs,wk(izt),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(izt),wk(iv1))
  call subqr (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(izt))
  go to 10
!
!  begin iteration loop.
!
! determine whether or not to stop --
!
 10   call inithv (1)
  nwpstp = nw - (iwfree-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,wk(izt),wk(iwfree),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iwfree-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
! form iterate.
!
  call vtriad (n,u,u,1e0/extrap,wk(izt),1)
!
! form residuals, as necessary.
!
  go to (161,162,163,164),iqlr + 1
!
 161  call suba (coef,jcoef,wfac,jwfac,n,u,wk(izt))
  call vexopy (n,wk(izt),rhs,wk(izt),2)
  go to 110
!
 162  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(izt))
  go to 110
!
 163  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subqr (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(izt))
  go to 110
!
 164  call suba (coef,jcoef,wfac,jwfac,n,u,wk(izt))
  call vexopy (n,wk(izt),rhs,wk(izt),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(izt),wk(iv1))
  call subqr (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(izt))
  go to 110
!
!  proceed to next iteration
!
 110  in = in + 1
  is = is + 1
  go to 10
!
!  finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'basicw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' basic method converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 996  call ershow (ier,'basicw')
  go to 735
!
! insuff. real wksp.
 999  ier = -2
  call ershow (ier,'basicw')
  go to 735
end
subroutine bbs (ndim,nn,maxt,t,x)
!
!*******************************************************************************
!
!! BBS does a banded back substitution.
!
!
!  (i + t)*x = y.
!
!
!     t is a rectangular matrix of adjacent super-diagonals.
!
!  Parameters:
!
!          ndim   row dimension of t array in defining routine
!          n      order of system
!          maxt   number of columns in t array
!          t      array of active size n by maxt giving the super-
!                  diagonals in the order 1,2,3,...
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension t(ndim,1), x(1)
!
  n = nn
  do 20 i = n-1,1,-1
     sum = x(i)
     lim = min (maxt,n-i)
     do 15 j = 1,lim
        sum = sum - t(i,j)*x(i+j)
 15      continue
     x(i) = sum
 20   continue
  return
end
subroutine bbsm (nsize,nsys,maxt,t,x)
!
!*******************************************************************************
!
!! BBSM does a banded back solve.
!
!
!  (i + t)*x = y.
!
!
!  T is an array containing superdiagonals in order 1,2,... .
!
!  Parameters:
!
!          n      order of system
!          nsize  size of a single subsystem
!          nsys   number of independent subsystems
!          maxt   number of columns in t array
!          t      array of active size n by maxt containing
!                  the super-diagonal elements of the factorization
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension t(nsize,nsys,1), x(nsize,1)
!
  do 25 i = nsize-1,1,-1
     lim = min (nsize-i, maxt)
     do 20 j = 1,lim
        ij = i + j
        do 15 l = 1,nsys
 15         x(i,l) = x(i,l) - t(i,l,j)*x(ij,l)
 20      continue
 25   continue
  return
end
subroutine bbst (ndim,nn,maxb,b,x)
!
!*******************************************************************************
!
!! BBST does a banded backward substitution.
!
!
!  (i + (b**t))*x = y.
!
!
!  The array b represents sub-diagonals.  b corresponds
!     to a banded system.
!
!  Parameters:
!
!          ndim   row dimension of b in defining routine
!          n      order of system (= nn)
!          maxb   number of diagonals stored in b
!          b      array of active size n x maxb giving the
!                  sub-diagonals in the order -1,-2,... .
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension b(ndim,1), x(1)
!
  n = nn
  do 25 i = n,2,-1
     term = x(i)
     lim = min (i-1,maxb)
     do 20 j = 1,lim
        x(i-j) = x(i-j) - b(i,j)*term
 20      continue
 25   continue
  return
end
subroutine bbstm (nsize,nsys,maxb,b,x)
!
!*******************************************************************************
!
!! BBSTM does the backward solve.
!
!
!  (i + (b**t))*x = y.
!
!
!  B contains subdiagonals for multiple banded systems.
!
!  Parameters:
!
!          n      order of system
!          nsize  the size of an individual subsystem
!          nsys   the number of subsystems
!          maxb   number of columns in b array
!          b      array of active size n by maxb containing
!                  sub-diagonals in the order -1,-2,-3,...
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension b(nsize,nsys,1), x(nsize,1)
!
  do 25 i = nsize,2,-1
     lim = min (i-1,maxb)
     do 20 j = 1,lim
        do 15 l = 1,nsys
 15         x(i-j,l) = x(i-j,l) - b(i,l,j)*x(i,l)
 20      continue
 25   continue
  return
end
subroutine bcgs (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef, &
  n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BCGS is the user interface to the biconjugate-gradient-squared algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2)
  dimension wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call bcgsw (suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs, &
    wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine bcgsw (suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
  wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BCGSW runs the biconjugate-gradient-squared algorithm.
!
!
! the algorithm is taken from "preconditioned biconjugate gradient
! methods for numerical reservoir simulation", by p. joly and r.
! eymard, to appear in journal of computational physics.  the original
! reference  is p. sonneveld, "cgs, a fast lanczos-type solver for
! unsymmetric linear systems," report 84-16, delft university of
! technology, dept. of mathematics and informatics.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subql, subqr
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, &
                    alphao, gamma, sigma, rr, rho, dkq, dkm1, &
                    ff, rqmin, rqmax, stptst, udnm, ubarnm, &
                    bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, &
         udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav, &
         rdot, rzdot, rztdot, zdot, zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
  nwusd = 0
  ier = 0
  iacel = 15
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  if (level >= 2) write (nout,496)
496   format (' bcgs')
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  if (iqr) go to 995
!
! initialize the stopping test.
!
  call inithv (0)
  zhave  = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
! allocate memory -- overlap wherever possible.
  ir0 = 1
  ip = ir0 + n
  ipt = ip + n
  if (.not. iqr) ipt = ip
  iq = ipt + n
  iz = iq + n
  izt = iz + n
  if (.not. iqr) izt = iz
  iv1 = izt + n
  iv2 = iv1 + n
  iv3 = iv2 + n
  nwusd = max(nwusd,iv3-1+n)
  ipaaq  = iv1
  ippaaq = iv2
!
! check the memory usage.
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  if (.not. iql) go to 121
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iz))
  go to 122
 121  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iz))
  call vexopy (n,wk(iz),rhs,wk(iz),2)
 122  if (iqr) call subqr (coef,jcoef,wfac,jwfac,n,wk(iz),wk(izt))
!
!  Begin iteration loop.
!
! determine whether or not to stop.
!
 10   call inithv (1)
  nwpstp = nw - (iv2-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,wk(iz),wk(izt),wk(iv2),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv2-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
  if (in /= 0) go to 110
!
! perform first-iterate calculations
!
  call vcopy (n,wk(iz),wk(ir0))
  call vcopy (n,wk(iz),wk(ip))
  call vcopy (n,wk(iz),wk(iq))
  r0r = vdot (n,wk(iz),wk(ir0))
  go to 111
!
! perform subsequent-iterate calculations
!
 110  r0rold = r0r
  r0r = vdot (n,wk(ir0),wk(iz))
  if (abs(r0rold) < srelpr**2) go to 996
  beta = r0r/r0rold
!
! form direction vectors.
!
  call vtriad (n,wk(ip),wk(iz),beta,wk(ipaaq),1)
  call vtriad (n,wk(iv2),wk(ipaaq),beta,wk(iq),1)
  call vtriad (n,wk(iq),wk(ip),beta,wk(iv2),1)
!
!  Form the iterate.
!
! at this point we have the vectors p and q and the new dot(r,r0).
! now form aq.
!
 111  iaq = iv1
  if (.not.iql) then
    call suba (coef,jcoef,wfac,jwfac,n,wk(iq),wk(iaq))
  else
    call suba (coef,jcoef,wfac,jwfac,n,wk(iq),wk(iv2))
    call subql (coef,jcoef,wfac,jwfac,n,wk(iv2),wk(iaq))
  end if
! dot(r0,aq).
  r0aq = vdot (n,wk(ir0),wk(iaq))
  if (abs(r0aq) < srelpr**2) go to 998
  alpha = r0r / r0aq
! p-alpha*aq, p+p-alpha*aq.
  call vtriad (n,wk(ipaaq), wk(ip),-alpha,wk(iaq),1)
  call vexopy (n,wk(ippaaq),wk(ip),wk(ipaaq),1)
!
!  get u.
  call vtriad (n,u,u,alpha,wk(ippaaq),1)
!
!  get resid.
  if (.not.iql) then
    call suba (coef,jcoef,wfac,jwfac,n,wk(ippaaq),wk(iv3))
    call vtriad (n,wk(iz),wk(iz),-alpha,wk(iv3),1)
  else
    call suba  (coef,jcoef,wfac,jwfac,n,wk(ippaaq),wk(iv3))
    call subql (coef,jcoef,wfac,jwfac,n,wk(iv3),wk(iv2))
    call vtriad (n,wk(iz),wk(iz),-alpha,wk(iv2),1)
  end if
!
!  proceed to next iteration
!
  in = in + 1
  is = is + 1
  go to 10
!
!  finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'bcgsw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' bcgs converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 995  ier = -16
  call ershow (ier,'bcgsw')
  go to 725
!
 996  ier = -13
  call ershow (ier,'bcgsw')
  go to 725
!
 997  call ershow (ier,'bcgsw')
  go to 735
!
 998  ier = -15
  call ershow (ier,'bcgsw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'bcgsw')
  go to 735
!
end
subroutine bdfac (lda,nn,nsizee,nt,nb,a,isym)
!
!*******************************************************************************
!
!! BDFAC computes the factorization of a dense banded matrix.
!
!
!  Parameters:
!
!        lda    leading dimension of array a
!        n      active size of array a
!        nsize  size of an individual subsystem (if multiple systems)
!                nsize = n upon input if not a multiple system
!        nt     number of diagonals needed to store the super-
!                diagonals
!        nb     number of diagonals needed to store the sub-
!                diagonals
!        a      array
!        isym   symmetry switch
!                = 0   matrix is symmetric
!                = 1   matrix is nonsymmetric
!
!  
!
  dimension a(lda,5)
  data lenv / 10 /
!
  n = nn
  maxt = nt
  nsize = nsizee
  nsys = n/nsize
!
!  branch on symmetry.
!
  if (isym == 1) go to 30
!
!  symmetric case.
!
!  diagonal case (maxt = 0).
!
  if (maxt /= 0) go to 15
  do 10 i = 1,n
 10   a(i,1) = 1.0/a(i,1)
  return
!
!  tridiagonal case (maxt = 1).
!
 15   if (maxt /= 1) go to 20
  if (nsys <= lenv) call tfac (n,a,a(1,2))
  if (nsys > lenv) call tfacm (n,nsize,a,a(1,2))
  return
!
!  pentadiagonal case (maxt = 2).
!
 20   if (maxt /= 2) go to 25
  if (nsys <= lenv) call pfac (n,a,a(1,2),a(1,3))
  if (nsys > lenv) call pfacm (n,nsize,a,a(1,2),a(1,3))
  return
!
!  banded case (maxt > 2).
!
 25   if (nsys <= lenv) call bfac (lda,n,maxt,a,a(1,2))
  if (nsys > lenv) call bfacm (n,nsize,nsys,maxt,a,a(1,2))
  return
!
!  nonsymmetric case.
!
 30   maxb = nb
!
!  diagonal case (maxt = maxb = 0).
!
  if (maxt /= 0 .or. maxb /= 0) go to 40
  do 35 i = 1,n
 35   a(i,1) = 1.0/a(i,1)
  return
!
!  tridiagonal case (maxt = maxb = 1).
!
 40   if (maxt /= 1 .or. maxb /= 1) go to 45
  if (nsys <= lenv) call tfacn (n,a,a(1,2),a(2,3))
  if (nsys > lenv) call tfacnm (n,nsize,a,a(1,2),a(2,3))
  return
!
!  pentadiagonal case (maxt = maxb = 2).
!
 45   if (maxt /= 2 .or. maxb /= 2) go to 50
  if (nsys <= lenv) call pfacn (n,a,a(1,2),a(1,3),a(2,4),a(3,5))
  if (nsys > lenv) call pfacnm (n,nsize,a,a(1,2),a(1,3),a(2,4),a(3,5))
  return
!
!  all other cases.
!
 50   if (nsys <= lenv) call bfacn (lda,n,maxt,maxb,a,a(1,2),a(1,maxt+2))
  if (nsys > lenv) call bfacnm (n,nsize,nsys,maxt,maxb,a,a(1,2),a(1,maxt+2))
  return
end
subroutine bdinv (lda,nn,nsizee,nt,nb,fac,isym)
!
!*******************************************************************************
!
!! BDINV computes the inverse of a dense banded matrix.
!
!
!  Parameters:
!
!        lda    leading dimension of factorization matrix fac
!        n      active size of factorization matrix fac
!        nsize  size of an individual subsystem (if multiple systems)
!                nsize = n upon input if not a multiple system
!        nt     number of diagonals needed to store the super-
!                diagonals
!        nb     number of diagonals needed to store the sub-
!                diagonals
!        fac    array containing factorization upon input
!        isym   symmetry switch
!                = 0   matrix is symmetric
!                = 1   matrix is nonsymmetric
!
!  
!
  dimension fac(lda,3)
  data lenv / 10 /
!
  n = nn
  maxt = nt
  nsize = nsizee
  nsys = n/nsize
!
!  branch on symmetry.
!
  if (isym == 1) go to 30
!
!  symmetric case.
!
  if (maxt - 1) 10,20,25
!
!  diagonal case (maxt = 0).
!
 10   return
!
!  tridiagonal case (maxt = 1).
!
 20   if (nsys <= lenv) call tinv (n,fac,fac(1,2))
  if (nsys > lenv) call tinvm (n,nsize,fac,fac(1,2))
  return
!
!  banded case (maxt >= 2).
!
 25   call binv (lda,n,maxt+1,fac)
  return
!
!  nonsymmetric case.
!
 30   maxb = nb
!
!  diagonal case (maxt = maxb = 0).
!
  if (maxt /= 0 .or. maxb /= 0) go to 40
  return
!
!  tridiagonal case (maxt = maxb = 1).
!
 40   if (maxt /= 1 .or. maxb /= 1) go to 45
  if (nsys <= lenv) call tinvn (n,fac,fac(1,2),fac(2,3))
  if (nsys > lenv) call tinvnm (n,nsize,fac,fac(1,2),fac(2,3))
  return
!
!  all other cases.
!
 45   call binvn (lda,n,maxt,maxb,fac,fac(1,2),fac(1,maxt+2))
  return
end
subroutine bdsol (lda,nn,nsizee,nt,nb,fac,y,x,isym)
!
!*******************************************************************************
!
!! BDSOL computes the solution to a dense banded matrix.
!
!
!     thus, bdsol finds the solution to   a*x = y,  where fac
!     contains the factorization of the a matrix.
!
!  Parameters:
!
!        lda    leading dimension of array fac
!        n      active size of array fac
!        nsize  size of an individual subsystem (if multiple systems)
!                nsize = n upon input if not a multiple system
!        nt     number of diagonals needed to store the super-
!                diagonals of the factorization
!        nb     number of diagonals needed to store the sub-
!                diagonals of the factorization
!        fac    array containing the factorization of the matrix
!        y      upon input, y conains the right hand side
!        x      upon output, x contains the solution to  a*x = y
!        isym   symmetry switch
!                = 0   matrix is symmetric
!                = 1   matrix is nonsymmetric
!
!  
!
  dimension fac(lda,5), x(1), y(1)
  data lenv / 10 /
!
  n = nn
  maxt = nt
  nsize = nsizee
  nsys = n/nsize
!
!  branch on symmetry.
!
  if (isym == 1) go to 30
!
!  symmetric case.
!
!  diagonal case (maxt = 0).
!
  if (maxt /= 0) go to 15
  do 10 i = 1,n
 10   x(i) = fac(i,1)*y(i)
  return
!
!  tridiagonal case (maxt = 1).
!
 15   if (maxt /= 1) go to 20
  if (nsys <= lenv) call tsoln (n,fac,fac(1,2),fac(1,2),y,x)
  if (nsys > lenv) call tsolnm (n,nsize,fac,fac(1,2),fac(1,2),y,x)
  return
!
!  pentadiagonal case (maxt = 2).
!
 20   if (maxt /= 2) go to 25
  if (nsys <= lenv) call psoln (n,fac,fac(1,2),fac(1,3),fac(1,2),fac(1,3),y,x)
  if (nsys > lenv) then
    call psolnm (n,nsize,fac,fac(1,2),fac(1,3),fac(1,2),fac(1,3),y,x)
  end if

  return
!
!  banded case (maxt >= 3).
!
 25   if (nsys <= lenv) call bsol (lda,n,maxt,fac,fac(1,2),y,x)
  if (nsys > lenv) call bsolm (n,nsize,maxt,fac,fac(1,2),y,x)
  return
!
!  nonsymmetric case.
!
 30   maxb = nb
!
!  diagonal case (maxt = maxb = 0).
!
  if (maxt /= 0 .or. maxb /= 0) go to 40
  do 35 i = 1,n
 35   x(i) = fac(i,1)*y(i)
  return
!
!  tridiagonal case (maxt = maxb = 1).
!
 40   if (maxt /= 1 .or. maxb /= 1) go to 45
  if (nsys <= lenv) call tsoln (n,fac,fac(1,2),fac(2,3),y,x)
  if (nsys > lenv) call tsolnm (n,nsize,fac,fac(1,2),fac(2,3), y,x)
  return
!
!  pentadiagonal case (maxt = maxb = 2).
!
 45   if (maxt /= 2 .or. maxb /= 2) go to 50
  if (nsys <= lenv) call psoln (n,fac,fac(1,2),fac(1,3),fac(2,4),fac(3,5),y,x)
  if (nsys > lenv) then
    call psolnm (n,nsize,fac,fac(1,2),fac(1,3),fac(2,4),fac(3,5),y,x)
  end if
  return
!
!  all other cases.
!
 50   if (nsys <= lenv) then
        call bsoln (lda,n,maxt,maxb,fac,fac(1,2),fac(1,maxt+2),y,x)
      end if

  if (nsys > lenv) then
    call bsolnm (n,nsize,maxt,maxb,fac,fac(1,2),fac(1,maxt+2),y,x)
  end if

  return
end
subroutine bdsolt (lda,nn,nsizee,nt,nb,fac,y,x)
!
!*******************************************************************************
!
!! BDSOLT computes the transpose solution to a nonsymmetric dense banded matrix.
!
!
!     thus, bdsolt finds the solution to   (a**t)*x = y,  where fac
!     contains the factorization of the a matrix.
!
!  Parameters:
!
!        lda    leading dimension of array fac
!        n      active size of array fac
!        nsize  size of an individual subsystem (if multiple systems)
!                nsize = n upon input if not a multiple system
!        nt     number of diagonals needed to store the super-
!                diagonals of the factorization
!        nb     number of diagonals needed to store the sub-
!                diagonals of the factorization
!        fac    array containing the factorization of the matrix
!        y      upon input, y conains the right hand side
!        x      upon output, x contains the solution to  a*x = y
!
!  
!
  dimension fac(lda,5), x(1), y(1)
  data lenv / 10 /
!
  n = nn
  maxt = nt
  maxb = nb
  nsize = nsizee
  nsys = n/nsize
!
!  nonsymmetric case.
!
!  diagonal case (maxt = maxb = 0).
!
  if (maxt /= 0 .or. maxb /= 0) go to 15
  do 10 i = 1,n
 10   x(i) = fac(i,1)*y(i)
  return
!
!  tridiagonal case (maxt = maxb = 1).
!
 15   if (maxt /= 1 .or. maxb /= 1) go to 20
  if (nsys <= lenv) call tsoln (n,fac,fac(2,3),fac(1,2),y,x)
  if (nsys > lenv) call tsolnm (n,nsize,fac,fac(2,3),fac(1,2),y,x)
  return
!
!  pentadiagonal case (maxt = maxb = 2).
!
 20   if (maxt /= 2 .or. maxb /= 2) go to 25
  if (nsys <= lenv) then
    call psoln (n,fac,fac(2,4),fac(3,5),fac(1,2),fac(1,3),y,x)
  end if

  if (nsys > lenv) then
    call psolnm (n,nsize,fac,fac(2,4),fac(3,5),fac(1,2),fac(1,3),y,x)
  end if

  return
!
!  all other cases.
!
 25   if (nsys <= lenv) then
    call bsolnt (lda,n,maxt,maxb,fac,fac(1,2),fac(1,maxt+2),y,x)
  end if

  if (nsys > lenv) then
    call bsontm (n,nsize,maxt,maxb,fac,fac(1,2),fac(1,maxt+2),y,x)
  end if

  return
end
subroutine bfac (ndim,nn,maxt,d,t)
!
!*******************************************************************************
!
!! BFAC computes a factorization to a single banded symmetric matrix.
!
!
!  Parameters:
!
!          ndim   row dimension of t array in defining routine
!          n      order of system (= nn)
!          maxt   number of columns in t array
!          d      vector containing the diagonal elements of a
!          t      array of active size n by maxt containing the
!                  super-diagonals in the order 1,2,3,...
!
  dimension d(1), t(ndim,1)
!
  n = nn
  nm1 = n - 1
  do 20 k = 1,nm1
     pivot = d(k)
     lim = min (n-k,maxt)
     do 15 j1 = 1,lim
        term = t(k,j1)/pivot
        jcol1 = k + j1
        d(jcol1) = d(jcol1) - term*t(k,j1)
        if (j1 == lim) go to 15
        j1p1 = j1 + 1
        do 10 j2 = j1p1,lim
           jcol2 = j2 - j1
           t(jcol1,jcol2) = t(jcol1,jcol2) - term*t(k,j2)
 10         continue
 15      continue
 20   continue
  do 25 i = 1,n
 25   d(i) = 1.0/d(i)
  do 35 j = 1,maxt
     len = n - j
     do 30 i = 1,len
 30      t(i,j) = d(i)*t(i,j)
 35   continue
  return
end
subroutine bfacm (n,nsize,nsys,maxt,d,t)
!
!*******************************************************************************
!
!! BFACM computes factorizations to multiple banded symmetric matrices.
!
!
!  Parameters:
!
!          n      order of global system (= nn)
!          nsize  order of a single system
!          nsys   number of independent subsystems
!          maxt   number of columns in t array
!          d      vector of length n containing the diagonal
!                  elements of a
!          t      array of active size n by maxt containing the
!                  super-diagonals in the order 1,2,3,...
!
!  
!
  dimension d(nsize,1), t(nsize,nsys,1)
!
  nsm1 = nsize - 1
  do 30 k = 1,nsm1
     lim = min (nsize-k,maxt)
     do 25 j1 = 1,lim
        jcol1 = k + j1
        do 10 l = 1,nsys
 10         d(jcol1,l) = d(jcol1,l) - (t(k,l,j1)**2)/d(k,l)
        if (j1 == lim) go to 25
        j1p1 = j1 + 1
        do 20 j2 = j1p1,lim
           jcol2 = j2 - j1
           do 15 l = 1,nsys
              t(jcol1,l,jcol2) = t(jcol1,l,jcol2) - t(k,l,j1)*t(k,l,j2)/d(k,l)
 15            continue
 20         continue
 25      continue
 30   continue
  call vinv (n,d)
  do 35 jj = 1,maxt
     len = n - jj
     call vexopy (len,t(1,1,jj),d,t(1,1,jj),3)
 35   continue
  return
end
subroutine bfacmy (methf,factor,coef,jcoef,wksp,iwksp,nn,ier)
!
!*******************************************************************************
!
!! BFACMY computes a block factorization.  (multicolor nonsymmetric diagonal)
!
!
!  parameters
!
!       n       order of system
!       nfactr  amount of real workspace needed for factorization
!       ier     error flag
!
!  
!
!
  external factor
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, &
                    ifact, kblsz, lvfill, ltrunc, ndeg, &
                    ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  dimension coef(1), wksp(1)
  integer jcoef(2), iwksp(1)
!
  n = nn
  if (methf <= 2) ivers = 1
  if (methf > 2) ivers = 2
!
!  calculate constants.
!
  if (ipropa == 0) propa = .false.
  if (ipropa == 1) propa = .true.
!
!  calculate fill-in and factor.
!
  call fillbc (n,ncolor,coef,jcoef,iwksp(iblock),wksp,iwksp,ier)
  if (ier < 0) return
  nwdiag = ndt + ndb + 1
  nwnew = nwdiag + 2*ltrunc
  if (methf == 1) nwkp = ncmax*nwnew
  if (methf == 2) nwkp = ncmax*(nwnew + 1)
  if (methf == 3) nwkp = 0
  if (methf == 4) nwkp = n + 2*ncmax
  call needw ('bfacmy',0,irpnt,nwkp,ier)
  if (ier < 0) return
  if (propa) then
     call factor (n,ndim,n,iwksp(iipnt),iwksp(jcnew+ncolor*nwdiag), &
       wksp(ifactr),coef(ndim*nwdiag+1),ncolor, &
       iwksp(nc),iwksp(iblock),iwksp(lbhb),0,1, &
       iwksp(ipt),omega,wksp(irpnt),ier)
  end if
  if (.not. propa) then
     call factor (n,n,n,iwksp(iipnt),iwksp(jcnew+ncolor*nwdiag), &
       wksp(ifactr),wksp(iwkpt2),ncolor,iwksp(nc),iwksp(iblock), &
       iwksp(lbhb),0,0,iwksp(ipt),omega,wksp(irpnt),ier)
  end if
  return
end
subroutine bfacmz (methf,factor,coef,jcoef,wksp,iwksp,nn,ier)
!
!*******************************************************************************
!
!! BFACMZ computes a block factorization.  (nonsymmetric diagonal)
!
!
!  parameters
!
!       n       order of system
!       nfactr  amount of real workspace needed for factorization
!       ier     error flag
!
!  
!
!
  external factor
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, &
         iblock, ncmax
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  dimension coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
!
  n = nn
  if (methf <= 2) ivers = 1
  if (methf > 2) ivers = 2
!
!  if requested, find out if matrix has block property a.
!
  ncol = n/kblsz
  if (ipropa == 0) propa = .false.
  if (ipropa == 1) propa = .true.
  if (lvfill > 0) propa = .false.
  if (lvfill > 0) go to 15
  if (ipropa /= 2) go to 15
  call needw ('bfacmz',1,iipnt,2*ncol,ier)
  if (ier < 0) return
  iwksp(iipnt) = lbhb
  call prbblk (ncol,1,iwksp(iblock),iwksp(iipnt),iwksp(iipnt+1), &
    iwksp(iipnt+ncol+1),propa)
  if (propa) ipropa = 1
  if (.not. propa) ipropa = 0
!
!  calculate fill-in and factor.
!
 15   call fillbn (n,coef,jcoef,iwksp(iblock),wksp,iwksp,ier)
  if (ier < 0) return
  nwnew = iwksp(iblock+2) + iwksp(iblock+5)
  nwdiag = nwnew - 2*ltrunc
  if (methf == 1) nwkp = kblsz*nwnew
  if (methf == 2) nwkp = kblsz*(nwnew + 1)
  if (methf == 3) nwkp = 0
  if (methf == 4) nwkp = n + 2*kblsz
  call needw ('fillbn',0,irpnt,nwkp,ier)
  if (ier < 0) return
  ipt1 = iblock + 3*lbhb
  ipt2 = ipt1 + nwnew
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
  if (propa) then
     call factor (n,ndim,n,iwksp(iipnt),jcoef(nwdiag+1),wksp(ifactr), &
       coef(ndim*nwdiag+1),1,idumb(1),iwksp(iblock),idumb(3),1,1, &
       idumb(2),omega,wksp(irpnt),ier)
  end if
  if (.not. propa .and. lvfill == 0) then
     call factor (n,n,n,iwksp(iipnt),jcoef(nwdiag+1),wksp(ifactr), &
       wksp(iwkpt2),1,idumb(1),iwksp(iblock),idumb(3),1,0, &
       idumb(2),omega,wksp(irpnt),ier)
  end if

  if (lvfill > 0) then
     call factor (n,n,n,iwksp(ipt1),iwksp(ipt2),wksp(ifactr),wksp(iwkpt2),1, &
       idumb(1),iwksp(iblock),idumb(3),1,0,idumb(2),omega,wksp(irpnt),ier)
  end if
  return
end
subroutine bfacn (ndim,nn,maxt,maxb,d,t,b)
!
!*******************************************************************************
!
!! BFACN computes a factorization to a single banded nonsymmetric matrix.
!
!
!  Parameters:
!
!          ndim   row dimension of t and b in defining routine
!          n      order of system (= nn)
!          maxt   number of diagonals stored in t
!          maxb   number of diagonals stored in b
!          d      vector of length n containing the diagonal
!                  elements of a
!          t      array of active size n x maxt giving the
!                  super-diagonals in the order 1,2,3,...
!          b      array of active size n x maxb giving the
!                  sub-diagonals in the order -1,-2,-3,...
!
!  
!
  dimension d(1), t(ndim,1), b(ndim,1)
!
  n = nn
  nm1 = n - 1
  do 35 k = 1,nm1
     pivot = d(k)
     liml = min (maxb,n-k)
     limu = min (maxt,n-k)
     do 30 ip = 1,liml
        i = k + ip
        term = b(i,ip)/pivot
        do 25 jp = 1,limu
           term1 = term*t(k,jp)
           l = jp - ip
           if (l) 10,15,20
 10            b(i,-l) = b(i,-l) - term1
           go to 25
 15            d(i) = d(i) - term1
           go to 25
 20            t(i,l) = t(i,l) - term1
 25         continue
 30      continue
 35   continue
!
  do 40 i = 1,n
 40   d(i) = 1.0/d(i)
  do 50 j = 1,maxt
     len = n - j
     do 45 i = 1,len
 45      t(i,j) = d(i)*t(i,j)
 50   continue
  do 60 j = 1,maxb
     len = n - j
     do 55 i = 1,len
 55      b(i+j,j) = d(i)*b(i+j,j)
 60   continue
  return
end
subroutine bfacnm (nn,nsize,nsys,maxt,maxb,d,t,b)
!
!*******************************************************************************
!
!! BFACNM computes a factorization to multiple banded nonsymmetric matrices.
!
!
!  Parameters:
!
!          nsize  size of a subsystem
!          nsys   number of independent subsystems
!          maxt   number of diagonals stored in t
!          maxb   number of diagonals stored in b
!          n      order of system (= nn)
!          d      vector of length n containing the diagonal
!                  elements of a
!          t      array of active size n x maxt giving the
!                  super-diagonals in the order 1,2,3,...
!          b      array of active size n x maxb giving the
!                  sub-diagonals in the order -1,-2,-3,...
!
!  
!
  dimension d(nsize,1), t(nsize,nsys,1), b(nsize,nsys,1)
!
  n = nn
  nsm1 = nsize - 1
  do 50 k = 1,nsm1
     liml = min (maxb,nsize-k)
     limu = min (maxt,nsize-k)
     do 45 ip = 1,liml
        i = k + ip
        do 40 jp = 1,limu
           l = jp - ip
           if (l) 10,20,30
 10            do 15 m = 1,nsys
 15            b(i,m,-l) = b(i,m,-l) - b(i,m,ip)*t(k,m,jp)/d(k,m)
           go to 40
 20            do 25 m = 1,nsys
 25            d(i,m) = d(i,m) - b(i,m,ip)*t(k,m,jp)/d(k,m)
           go to 40
 30            do 35 m = 1,nsys
 35            t(i,m,l) = t(i,m,l) - b(i,m,ip)*t(k,m,jp)/d(k,m)
 40         continue
 45      continue
 50   continue
!
  call vinv (n,d)
  do 55 j = 1,maxt
     len = n - j
     call vexopy (len,t(1,1,j),d,t(1,1,j),3)
 55   continue
  do 60 j = 1,maxb
     len = n - j
     call vexopy (len,b(j+1,1,j),d,b(j+1,1,j),3)
 60   continue
  return
end
subroutine bfacs (methf,factor,coef,jcoef,wksp,iwksp,nn,ier)
!
!*******************************************************************************
!
!! BFACS computes a block factorization.  (symmetric diagonal)
!
!
!  parameters
!
!       n       order of system
!       nfactr  amount of real workspace needed for factorization
!       ier     error flag
!
!  
!
!
  external factor
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, &
         iblock, ncmax
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  integer jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
!
  n = nn
  if (methf <= 2) ivers = 1
  if (methf > 2) ivers = 2
!
!  if requested, find out if matrix has block property a.
!
  ncol = n/kblsz
  if (ipropa == 0) propa = .false.
  if (ipropa == 1) propa = .true.
  if (lvfill > 0) propa = .false.
  if (lvfill > 0) go to 15
  if (ipropa /= 2) go to 15
  call needw ('bfacs',1,iipnt,2*ncol,ier)
  if (ier < 0) return
  iwksp(iipnt) = lbhb
  call prbblk (ncol,1,iwksp(iblock),iwksp(iipnt), iwksp(iipnt+1), &
    iwksp(iipnt+ncol+1),propa)
  if (propa) ipropa = 1
  if (.not. propa) ipropa = 0
!
!  calculate fill-in and factor.
!
 15   call fillb (n,coef,jcoef,iwksp(iblock),wksp,iwksp,ier)
  if (ier < 0) return
  nwnew = iwksp(iblock+2)
  nwdiag = nwnew - ltrunc
  if (methf == 1) nwkp = kblsz*nwnew
  if (methf == 2) nwkp = kblsz*(nwnew + 1)
  if (methf == 3) nwkp = 0
  if (methf == 4) nwkp = n + 2*kblsz
  call needw ('fillb',0,irpnt,nwkp,ier)
  if (ier < 0) return
  ipt1 = iblock + 3*lbhb
  ipt2 = ipt1 + nwnew
  if (propa) then
     call factor (n,ndim,n,iwksp(iipnt),jcoef(nwdiag+1),wksp(ifactr), &
       coef(ndim*nwdiag+1),kblsz,iwksp(iblock),lbhb,1,omega,wksp(irpnt),ier)
  end if
  if (.not. propa .and. lvfill == 0) then
     call factor (n,n,n,iwksp(iipnt),jcoef(nwdiag+1),wksp(ifactr), &
       wksp(iwkpt2),kblsz,iwksp(iblock),lbhb,0,omega,wksp(irpnt),ier)
  end if
  if (lvfill > 0) then
     call factor (n,n,n,iwksp(ipt1),iwksp(ipt2),wksp(ifactr),wksp(iwkpt2), &
      kblsz,iwksp(iblock),lbhb,0,omega,wksp(irpnt),ier)
  end if
  return
end
subroutine bfs (ndim,nn,maxb,b,x)
!
!*******************************************************************************
!
!! BFS does a forward substitution.
!
!
!  (i + b)*x = y.
!
!
!  The array b represents sub-diagonals.  b corresponds to a banded system.
!
!  Parameters:
!
!          ndim   row dimension of b in defining routine
!          n      order of system (= nn)
!          maxb   number of diagonals stored in b
!          b      array of active size n x maxb giving the
!                  sub-diagonals in the order -1,-2,-3,... .
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension b(ndim,1), x(1)
!
  n = nn
  do 15 i = 2,n
     lim = min (i-1,maxb)
     sum = x(i)
     do 10 j = 1,lim
        sum = sum - b(i,j)*x(i-j)
 10      continue
     x(i) = sum
 15   continue
  return
end
subroutine bfsm (nsize,nsys,maxb,b,x)
!
!*******************************************************************************
!
!! BFSM does the forward solve.
!
!
!  (i + b)*x = y.
!
!
!  B contains subdiagonals for multiple banded systems.
!
!  Parameters:
!
!          n      order of system
!          nsize  the size of an individual subsystem
!          nsys   the number of subsystems
!          maxb   number of columns in b array
!          b      array of active size n by maxb containing
!                  sub-diagonals in the order -1,-2,-3,... .
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension b(nsize,nsys,1), x(nsize,1)
!
  do 20 i = 2,nsize
     lim = min (i-1,maxb)
     do 15 j = 1,lim
        do 10 l = 1,nsys
 10         x(i,l) = x(i,l) - b(i,l,j)*x(i-j,l)
 15      continue
 20   continue
  return
end
subroutine bfst (ndim,nn,maxt,t,x)
!
!*******************************************************************************
!
!! BFST does a banded forward substitution.
!
!
!   (i + (t**t))*x = y.
!
!
!     t is a rectangular matrix of adjacent super-diagonals.
!
!  Parameters:
!
!          ndim   row dimension of t array in defining routine
!          n      order of system
!          maxt   number of columns in t array
!          t      array of active size n by maxt giving the super-
!                  diagonals in the order 1,2,3,...
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension t(ndim,1), x(1)
!
  n = nn
  nm1 = n - 1
  do 20 i = 1,nm1
     term = x(i)
     lim = min (maxt,n-i)
     do 15 j = 1,lim
        x(i+j) = x(i+j) - t(i,j)*term
 15      continue
 20   continue
  return
end
subroutine bfstm (nsize,nsys,maxt,t,x)
!
!*******************************************************************************
!
!! BFSTM does a forward solve.
!
!
!  (i + (t**t))*x = y.
!
!
!  T is an array containing superdiagonals in order 1,2,... .
!     (multiple systems)
!
!  Parameters:
!
!          n      order of system
!          nsize  size of a single subsystem
!          nsys   number of independent subsystems
!          maxt   number of columns in t array
!          t      array of active size n by maxt containing
!                  the super-diagonal elements of the factorization
!          x      on input, x contains y
!                 vector containing solution upon output
!
!  
!
  dimension t(nsize,nsys,1), x(nsize,1)
!
  nsm1 = nsize - 1
  do 20 i = 1,nsm1
     lim = min (maxt,nsize-i)
     do 15 j = 1,lim
        ij = i + j
        do 10 l = 1,nsys
 10         x(ij,l) = x(ij,l) - t(i,l,j)*x(i,l)
 15      continue
 20   continue
  return
end
subroutine bic2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BIC2 drives the block factorization (version 1) method.
!
  external accel, suba1, subq25, copy, noadp
  external ibfcs1
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacs (1,ibfcs1,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  iwkpt1 = irpnt
  irpnt = irpnt + kblsz
  if (ier < 0) return
  call split (accel,suba1,suba1,subq25,subq25,subq25,subq25,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - kblsz
  return
end
subroutine bic3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BIC3 drives the block factorization (version 1) method.
!
  external accel, suba4, suba5, subq70, subq71, subq72, subq73, subq74
  external subq75, noadp
  external ibfcn1
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacmz (1,ibfcn1,coef,jcoef,wksp,iwksp, n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*kblsz
  call split (accel,suba4,suba5,subq70,subq71,subq72,subq73,subq74,subq75, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*kblsz
  return
end
subroutine bic7 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BIC7 drives the block factorization (version 1) method.  (multi-color ordering)
!
  external accel, suba2, suba3, subq34, subq35, subq36
  external subq37, subq38, subq39, noadp
  external ibfcn1
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, &
         iblock, ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  t1 = timer (dummy)
  if (ifact == 1) call bfacmy (1,ibfcn1,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*ncmax
  call split (accel,suba2,suba3,subq34,subq35,subq36,subq37,subq38,subq39, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*ncmax
  return
end
subroutine bicol (n,nz,ia,ja,count,father,oppos,propa)
!
!*******************************************************************************
!
!! BICOLOR determines whether or not a matrix is bi-colorable.
!
!
!  The matrix is represented in the sparse (ia,ja) format is bi-colorable.
!     the algorithm used is the union-find algorithm.
!
!  Parameters:
!
!        n      number of vertices
!        nz     number of edges (length of ia and ja vectors)
!        ia     integer vector of i values
!        ja     integer vector of j values
!        count  integer workspace vectors of length n each
!        father upon output, count gives the color of each node
!        oppos
!        propa  logical variable indicating on output whether
!                matrix has property a
!
!  specification of parameters
!
  logical propa
  integer ia(1), ja(1), count(1), father(1), oppos(1)
  integer v, w, w0, a, b, c, d
!
  do 10 i = 1,n
     count(i) = 1
     father(i) = 0
     oppos(i) = 0
 10   continue
  do 60 k = 1,nz
     if (ia(k) == ja(k)) go to 60
!
!  a = find (ia(k)).
!
     v = ia(k)
 15      if (father(v) == 0) go to 20
     v = father(v)
     go to 15
 20      w = ia(k)
 25      if (father(w) == 0) go to 30
     w0 = w
     w = father(w)
     father(w0) = v
     go to 25
 30      a = v
!
!  b = find (ja(k)).
!
     v = ja(k)
 35      if (father(v) == 0) go to 40
     v = father(v)
     go to 35
 40      w = ja(k)
 45      if (father(w) == 0) go to 50
     w0 = w
     w = father(w)
     father(w0) = v
     go to 45
 50      b = v
!
!  test for a = b.
!
     if (a /= b) go to 55
     propa = .false.
     return
!
!  do unioning.
!
 55      if (oppos(a) == b) go to 60
     if (oppos(b) == 0) then
        c = a
     else
!
!  c = merge (a,oppos(b)).
!
        i = a
        j = oppos(b)
        if (count(i) >= count(j)) then
           father(j) = i
           count(i) = count(i) + count(j)
           c = i
        else
           father(i) = j
           count(j) = count(i) + count(j)
           c = j
        end if
     end if
     if (oppos(a) == 0) then
        d = b
     else
!
!  d = merge (b,oppos(a)).
!
        i = b
        j = oppos(a)
        if (count(i) >= count(j)) then
           father(j) = i
           count(i) = count(i) + count(j)
           d = i
        else
           father(i) = j
           count(j) = count(i) + count(j)
           d = j
        end if
     end if
     oppos(c) = d
     oppos(d) = c
 60   continue
!
!  do coloring.
!
  do 65 i = 1,n
 65   count(i) = 0
  do 90 i = 1,n
!
!  a = find(i).
!
     v = i
 70      if (father(v) == 0) go to 75
     v = father(v)
     go to 70
 75      w = i
 80      if (father(w) == 0) go to 85
     w0 = w
     w = father(w)
     father(w0) = v
     go to 80
 85      a = v
     if (count(a) == 0) then
        count(a) = 1
        count(i) = 1
        j = oppos(a)
        if (j /= 0) count(j) = 2
     else
        count(i) = count(a)
     end if
 90   continue
  propa = .true.
  return
end
subroutine bicx2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BICX2 drives the block factorization (version 2) method.
!
  external accel, suba1, subq25, copy, noadp
  external ibfcs2
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacs (3,ibfcs2,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  iwkpt1 = irpnt
  irpnt = irpnt + kblsz
  if (ier < 0) return
  call split (accel,suba1,suba1,subq25,subq25,subq25,subq25,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - kblsz
  return
end
subroutine bicx3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BICX3 drives the block factorization (version 2) method.
!
  external accel, suba4, suba5, subq70, subq71, subq72
  external subq73, subq74, subq75, noadp
  external ibfcn2
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacmz (3,ibfcn2,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*kblsz
  call split (accel,suba4,suba5,subq70,subq71,subq72,subq73,subq74,subq75, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*kblsz
  return
end
subroutine bicx7 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! BICX7 drives the block factorization (version 2) method (multi-color ordering)
!
  external accel, suba2, suba3, subq34, subq35, subq36
  external subq37, subq38, subq39, noadp
  external ibfcn2
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, &
         iblock, ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  t1 = timer (dummy)
  if (ifact == 1) call bfacmy (3,ibfcn2,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*ncmax
  call split (accel,suba2,suba3,subq34,subq35,subq36,subq37,subq38,subq39, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*ncmax
  return
end
subroutine binv (ndim,nn,maxnz,fact)
!
!*******************************************************************************
!
!! BINV computes an approximate inverse to a single banded symmetric matrix.  
!
!
!  fact must contain upon input the output from a factorization routine.
!
!  Parameters:
!
!          ndim   row dimension of fact in the defining routine
!          n      order of system (= nn)
!          maxnz  bandwidth of the factorization and inverse
!          fact   array containing factorization diagonals
!                  in the order 0,1,2,3,...
!
!  
!
  dimension fact(ndim,2)
!
  n = nn
  nm1 = n - 1
!
!  general banded matrix.
!
  do 25 ik = 1,nm1
     k = n - ik
     lim = min (ik+1,maxnz)
     sum1= 0.0
     do 15 i = 2,lim
        t1 = fact(k,i)
        sum2= 0.0
        do 10 j = 2,lim
           m1 = min (i,j)
           m2 = max (i,j)
           l1 = k + m1 - 1
           l2 = m2 - m1 + 1
           sum2 = sum2 - fact(k,j)*fact(l1,l2)
 10         continue
        fact(n,i) = sum2
        sum1 = sum1 - t1*sum2
 15      continue
     fact(k,1) = fact(k,1) + sum1
     do 20 i = 2,lim
 20      fact(k,i) = fact(n,i)
 25   continue
  do 30 i = 2,maxnz
 30   fact(n,i)= 0.0
  return
end
subroutine binvn (ndim,nn,maxt,maxb,d,t,b)
!
!*******************************************************************************
!
!! BINVN computes an approximate inverse to a single banded nonsymmetric matrix.  
!
!
!  d, t, and b must contain upon input
!     the output from a factorization routine.
!
!  Parameters:
!
!          ndim   row dimension of t and b in the defining routine
!          n      order of system (= nn)
!          maxt   number of columns in t
!          maxb   number of columns in b
!          d      vector of length n containing the diagonal
!                  elements of the factorization
!          t      array of active size n by maxt containing
!                  the superdiagonals of the factorization
!                  in the order 1,2,3,...
!          b      array of active size n by maxb containing
!                  the subdiagonals of the factorization
!                   in the order -1,-2,-3,....
!
!  
!
  dimension d(1), t(ndim,1), b(ndim,1)
!
  n = nn
  nm1 = n - 1
!
!  general banded matrix.
!
  do 75 ik = 1,nm1
     k = n - ik
!
!  copy kth row and column into wksp.
!
     limr = min (maxt,ik)
     limc = min (maxb,ik)
     do 10 j = 1,limr
 10      t(n,j) = t(k,j)
     do 15 j = 1,limc
 15      b(1,j) = b(k+j,j)
!
!  do computations for kth row.
!
     do 40 j = 1,limr
        sum= 0.0
        lim = min (limr,limc+j)
        do 35 i = 1,lim
           kpi = k + i
           l = i - j
           if (l) 20,25,30
 20            sum = sum - t(n,i)*t(kpi,-l)
           go to 35
 25            sum = sum - t(n,i)*d(kpi)
           go to 35
 30            sum = sum - t(n,i)*b(kpi,l)
 35         continue
        t(k,j) = sum
 40      continue
!
!  do computations for kth column.
!
     do 65 j = 1,limc
        sum= 0.0
        lim = min (limc,limr+j)
        kpj = k + j
        do 60 i = 1,lim
           kpi = k + i
           l = i - j
           if (l) 45,50,55
 45            sum = sum - b(1,i)*b(kpj,-l)
           go to 60
 50            sum = sum - b(1,i)*d(kpi)
           go to 60
 55            sum = sum - b(1,i)*t(kpj,l)
 60         continue
        b(kpj,j) = sum
 65      continue
!
!  compute kth diagonal element.
!
     sum = d(k)
     lim = min (limr,limc)
     do 70 j = 1,lim
 70      sum = sum - t(n,j)*b(k+j,j)
     d(k) = sum
 75   continue
!
!  zero out workspace rows.
!
  do 80 j = 1,maxt
 80   t(n,j)= 0.0
  do 85 j = 1,maxb
 85   b(1,j)= 0.0
  return
end
subroutine blkdef (coef,jcoef,wksp,iwksp,nn,ier)
!
!*******************************************************************************
!
!! BLKDEF defines various block constants for a constant block size matrix.
!
!
!  Parameters:
!
!        n        problem size
!
!  common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, &
         iblock, ncmax
  integer jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
!
  n = nn
!
  call needw ('blkdef',1,iipnt,3*(maxnz+1),ier)
  if (ier < 0) return
  call move5 (ndim,n,maxnz,jcoef,coef)
  if (ifact == 0) return
  ifacti = iipnt
  iblock = ifacti
  call defcon (ndim,n,maxnz,jcoef,coef,kblsz,iwksp(ifacti),lbhb)
  nfacti = 3*lbhb
  iipnt = ifacti + 3*lbhb
  return
end
subroutine bmul (ndim,n,maxt,d,t,x,y)
!
!*******************************************************************************
!
!! BMUL computes y = A*x, where A is a banded symmetric matrix.
!
!
!  Parameters:
!
!         ndim          row dimension of array t
!         n             order of matrix
!         maxt          number of columns in t
!         d             vector of length n giving the
!                        diagonal elements of a
!         t             array of size n by maxt giving the
!                        superdiagonals of a in the order
!                        1,2,....
!         x,y           vectors of order n
!
!  
!
  dimension d(1), t(ndim,1), x(1), y(1)
!
  do 10 i = 1,n
 10   y(i) = d(i)*x(i)
  if (maxt <= 0) return
  do 25 la = 1,maxt
     len = n - la
     do 15 i = 1,len
 15      y(i) = y(i) + t(i,la)*x(i+la)
     do 20 i = 1,len
 20      y(i+la) = y(i+la) + t(i,la)*x(i)
 25   continue
  return
end
subroutine bmuln (ndim,n,maxt,maxb,d,t,b,x,y)
!
!*******************************************************************************
!
!! BMULN computes y = A*x, where A is in nonsymmetric band format.
!
!
!  A is represented by arrays d, t, and b.
!
!  Parameters:
!
!         ndim          row dimension of arrays t and b
!         n             order of array a
!         maxt          number of columns in t array
!         maxb          number of columns in b array
!         d             vector of length n giving the diagonal
!                        elements of a
!         t             array of active size n by maxt giving
!                        the super-diagonals of a in the order
!                        1,2,3,...
!         b             array of active size n by maxb giving
!                        the sub-diagonals of a in the order
!                        -1,-2,-3,....
!         x,y           vectors of order n
!
!  
!
  dimension d(1), t(ndim,1), b(ndim,1), x(1), y(1)
!
  do 10 i = 1,n
 10   y(i) = d(i)*x(i)
  if (maxt < 1) go to 25
  do 20 j = 1,maxt
     len = n - j
     do 15 i = 1,len
 15      y(i) = y(i) + t(i,j)*x(i+j)
 20   continue
 25   if (maxb < 1) return
  do 35 j = 1,maxb
     len = n - j
     do 30 i = 1,len
 30      y(i+j) = y(i+j) + b(i+j,j)*x(i)
 35   continue
  return
end
subroutine bmulnt (ndim,n,maxt,maxb,d,t,b,x,y)
!
!*******************************************************************************
!
!! BMULNT computes y = (A**t)*x, where A is in nonsymmetric band format.
!
!
!  A is represented by d, t, and b.
!
!  Parameters:
!
!         ndim          row dimension of arrays t and b
!         n             order of array a
!         maxt          number of columns in t array
!         maxb          number of columns in b array
!         d             vector of length n giving the diagonal
!                        elements of a
!         t             array of active size n by maxt giving
!                        the super-diagonals of a in the order
!                        1,2,3,...
!         b             array of active size n by maxb giving
!                        the sub-diagonals of a in the order
!                        -1,-2,-3,...
!         x,y           vectors of order n
!
!  
!
  dimension d(1), t(ndim,1), b(ndim,1), x(1), y(1)
!
  do 10 i = 1,n
 10   y(i) = d(i)*x(i)
  if (maxt < 1) go to 25
  do 20 j = 1,maxt
     len = n - j
     do 15 i = 1,len
 15      y(i+j) = y(i+j) + t(i,j)*x(i)
 20   continue
 25   if (maxb < 1) return
  do 35 j = 1,maxb
     len = n - j
     do 30 i = 1,len
 30      y(i) = y(i) + b(i+j,j)*x(i+j)
 35   continue
  return
end
subroutine bsol (ndim,nn,maxt,d,t,y,x)
!
!*******************************************************************************
!
!! BSOL solves A*x = y for a banded and symmetric matrix A. 
!
!
!  d and t must contain upon input the factorization arrays from bfac.
!
!  Parameters:
!
!          ndim   row dimension of t array in defining routine
!          n      order of system
!          maxt   number of columns in t array
!          d      vector of length n containing the diagonal
!                  pivots of the factorization
!          t      array of active size n by maxt giving the super-
!                  diagonals of the factorization in the order
!                  1,2,3,...
!          y      right-hand-side vector
!          x      vector containing solution upon output
!
!  
!
  dimension t(ndim,1), x(1), y(1), d(1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call bfst (ndim,n,maxt,t,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call bbs (ndim,n,maxt,t,x)
  return
end
subroutine bsolm (nn,nsize,maxt,d,t,y,x)
!
!*******************************************************************************
!
!! BSOLM solves the system A*x = y where A is multiple symmetric banded matrices 
!
!
!  The factorizations are contained in d and t.
!
!  Parameters:
!
!          n      order of system
!          nsize  size of a single subsystem
!          maxt   number of columns in t array
!          d      vector of length n containing the diagonal
!                  elements of the factorization
!          t      array of active size n by maxt containing
!                  the super-diagonal elements of the factorization
!                  in the order 1,2,3,...
!          y      right-hand-side vector
!          x      vector containing solution upon output
!
!  
!
  dimension d(1), t(1), y(1), x(1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  nsys = n/nsize
  call bfstm (nsize,nsys,maxt,t,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call bbsm (nsize,nsys,maxt,t,x)
  return
end
subroutine bsoln (ndim,nn,maxt,maxb,d,t,b,y,x)
!
!*******************************************************************************
!
!! BSOLN solves A*x = y for a banded and nonsymmetric matrix.
!
!
!     d, t, and b must contain upon input the factorization arrays
!     from bfacn.
!
!  Parameters:
!
!          ndim   row dimension of t array in defining routine
!          n      order of system
!          maxt   number of columns in t array
!          maxb   number of columns in b array
!          d      vector of length n containing the diagonal
!                  pivots of the factorization
!          t      array of active size n by maxt giving the super-
!                  diagonals of the factorization in the order
!                  1,2,3,...
!          b      array of active size n by maxb giving the sub-
!                  diagonals of the factorization in the order
!                  -1,-2,-3,...
!          y      right-hand-side vector
!          x      vector containing solution upon output
!
!  
!
  dimension t(ndim,1), x(1), y(1), d(1), b(ndim,1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call bfs (ndim,n,maxb,b,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call bbs (ndim,n,maxt,t,x)
  return
end
subroutine bsolnm (nn,nsize,maxt,maxb,d,t,b,y,x)
!
!*******************************************************************************
!
!! BSOLNM solves A*x = y for a banded and nonsymmetric matrix.
!
!
!     d, t, and b must contain upon input the factorization arrays
!     from bfacnm.  (multiple systems)
!
!  Parameters:
!
!          n      order of system
!          nsize  size of an individual subsystem
!          maxt   number of columns in t array
!          maxb   number of columns in b array
!          d      vector of length n containing the diagonal
!                  pivots of the factorization
!          t      array of active size n by maxt giving the super-
!                  diagonals of the factorization in the order
!                  1,2,3,...
!          b      array of active size n by maxb giving the sub-
!                  diagonals of the factorization in the order
!                  -1,-2,-3,...
!          y      right-hand-side vector
!          x      vector containing solution upon output
!
!  
!
  dimension t(1), x(1), y(1), d(1), b(1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  nsys = n/nsize
  call bfsm (nsize,nsys,maxb,b,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call bbsm (nsize,nsys,maxt,t,x)
  return
end
subroutine bsolnt (ndim,nn,maxt,maxb,d,t,b,y,x)
!
!*******************************************************************************
!
!! BSOLNT solves (A**t)*x = y for a banded and nonsymmetric matrix.
! 
!
!  d, t, and b must contain upon input the
!     factorization arrays from bfacn.
!
!  Parameters:
!
!          ndim   row dimension of t array in defining routine
!          n      order of system
!          maxt   number of columns in t array
!          maxb   number of columns in b array
!          d      vector of length n containing the diagonal
!                  pivots of the factorization
!          t      array of active size n by maxt giving the super-
!                  diagonals of the factorization in the order
!                  1,2,3,...
!          b      array of active size n by maxb giving the sub-
!                  diagonals of the factorization in the order
!                  -1,-2,-3,...
!          y      right-hand-side vector
!          x      vector containing solution upon output
!
!  
!
  dimension t(ndim,1), x(1), y(1), d(1), b(ndim,1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call bfst (ndim,n,maxt,t,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call bbst (ndim,n,maxb,b,x)
  return
end
subroutine bsontm (nn,nsize,maxt,maxb,d,t,b,y,x)
!
!*******************************************************************************
!
!! BSONTM solves (A**t)*x = y for a banded and nonsymmetric matrix.  
!
!  d, t, and b must contain upon input the
!     factorization arrays from bfacnm.  (multiple systems)
!
!  Parameters:
!
!          n      order of system
!          nsize  size of an individual subsystem
!          maxt   number of columns in t array
!          maxb   number of columns in b array
!          d      vector of length n containing the diagonal
!                  pivots of the factorization
!          t      array of active size n by maxt giving the super-
!                  diagonals of the factorization in the order
!                  1,2,3,...
!          b      array of active size n by maxb giving the sub-
!                  diagonals of the factorization in the order
!                  -1,-2,-3,...
!          y      right-hand-side vector
!          x      vector containing solution upon output
!
!  
!
  dimension t(1), x(1), y(1), d(1), b(1)
!
  n = nn

  x(1:n) = y(1:n)
  nsys = n/nsize
  call bfstm (nsize,nsys,maxt,t,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call bbstm (nsize,nsys,maxb,b,x)
  return
end
subroutine cg (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n, &
  u,ubar,rhs,wksp,iwksp, iparm,rparm,ier)
!
!*******************************************************************************
!
!! CG is the user interface to the conjugate gradient algorithm.
!
  external suba, subat, subql, subqlt, subqr, subqrt, subadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, &
         iplr, iqlr, ntest, is, iacel, idgts, nbl1d, nbl2d
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  ier = 0
  call needw ( 'CG', 0, irpnt, 3*n+2*itmax, ier )

  if ( ier < 0 ) then
    return
  end if

  nw = lenr - irpnt + 1

  call cgw ( suba, subql, coef, jcoef, wksp, iwksp, n, u, ubar, rhs, &
    wksp(irpnt), nw, iparm, rparm, ier )

  irmax = irpnt + nw - 1

  return
end
subroutine cgcr (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef, &
  n,u,ubar,rhs,wk,iwk,iparm,rparm,ier)
!
!*******************************************************************************
!
!! CGCR implements the constrained residual method.
!
!
!  The CGCR method of j. r. wallis is coupled with truncated/restarted 
!  orthomin.  for further information about the algorithm, see 
!  "constrained residual acceleration of conjugate residual methods", 
!  by j. r. wallis,
! r. p. kendall and t. e. little of j. s. nolen and assocs. inc.;
! report spe 13536, society of petroleum engineers, 1985.
!
! right preconditioning only is allowed in this algorithm.
!
! unfortunately, this routine is limited -- all blocks must be the
! same size.  but the idea can be easily generalized.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wk(1), iwk(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
  external nullpl, cgcrpr
  logical ipl, ipr, iql, iqr
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, &
         iplr, iqlr, ntest, is, iacel, idgts, nbl1d, nbl2d
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!  data common blocks
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / ccgcr  / nblk, nband, ictac, ieta, ivcgcr
!
! time to proceed.
!
  if (nstore/=2 .and. nstore/=3) go to 998
!
  irpsav = irpnt
  iql = iqlr==1 .or. iqlr==3
  iqr = iqlr==2 .or. iqlr==3
  if (iql) go to 998
!
  ipl = .false.
  ipr = .true.
  iplr = 0
  if (ipl) iplr = iplr + 1
  if (ipr) iplr = iplr + 2
!
! form the c**(t)*a*c matrix
!
 1    if (nbl1d<=0 .or. nbl2d<=0) go to 998
  nbl0d = 1
  if (mod(nbl2d,nbl1d)/=0 .or. mod(nbl1d,nbl0d)/=0) go to 998
  nblk = n / nbl2d
  if (nblk == 1) nblk = n / nbl1d
  ictac = irpnt
  nwgb = lenr - ictac + 1
  ierpp = 0
  call getblk (coef,jcoef,n,nblk,nband,wk(ictac),nwgb,ierpp)
  irmax = max (irmax,ictac-1+nwgb)
  if (ierpp < 0) go to 999
  irpnt = ictac + nblk*nband
!
! perform first-iterate calculations
!
  ieta = irpnt
  ivcgcr = ieta + nblk
  iv2 = ivcgcr + n
  irmax = max(irmax,iv2-1+n)
  if (irmax > lenr) go to 997
!
  call suba (coef,jcoef,wk,iwk,n,u,wk(ivcgcr))
  call vexopy (n,wk(ivcgcr),rhs,wk(ivcgcr),2)
  call tmult (n,nblk,nband,wk(ictac),wk(ieta),wk(ivcgcr),wk(ivcgcr))
  call vexopy (n,u,u,wk(ivcgcr),1)
!
! pass it on to orthomin.
!
  irpnt = iv2
  nw = lenr - irpnt + 1
  call omingw (suba,subql,subqr,nullpl,cgcrpr,coef,jcoef,wk,iwk,n,u,ubar, &
    rhs,wk(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
!
  irpnt = irpsav
  return
!
! error returns.
!
! insuff. real workspace.
 997  ier = -2
  call ershow (ier,'cgcr')
  return
!
! unimplemented option.
 998  ier = -16
  call ershow (ier,'cgcr')
  return
!
! generic handler.
 999  ier = ierpp
  return
end
subroutine cgcrpr (coef,jcoef,wk,iwk,n,subql,suba,subqr,u,v)
!
!*******************************************************************************
!
!! CGCRPR is a right preconditioner routine to use with the CGCR method.
!
  dimension u(1), v(1), coef(1), jcoef(2), wk(1), iwk(1)
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, &
         iplr, iqlr, ntest, is, iacel, idgts, nbl1d, nbl2d
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / ccgcr  / nblk, nband, ictac, ieta, ivcgcr
!
  external subql, suba, subqr
!
! could bypass next line if subqr is just a copy.
  call subqr (coef,jcoef,wk,iwk,n,u,v)
  call suba (coef,jcoef,wk,iwk,n,v,wk(ivcgcr))
  call tmult (n,nblk,nband,wk(ictac),wk(ieta),wk(ivcgcr), wk(ivcgcr))
  call vexopy (n,v,v,wk(ivcgcr),2)
!
  return
end
subroutine cgnr (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n, &
  u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! CGNR is the user interfact to the conjugate gradient algorithm on the normal equations.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call cgnrw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wksp,iwksp,n, &
    u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine cgnrw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! CGNRW runs the conjugate gradient algorithm on the normal equations.
!
!
! in this variant, the residual of the original system is minimized
! per iteration.  currently, only left preconditioning is implemented.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subat, subql, subqlt, subqr, subqrt
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, &
         iplr, iqlr, ntest, is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, sigma, &
         rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, bnorm, &
         bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 5
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  if (iqr) go to 995
  if (level >= 2) write (nout,496)
496   format (' cgnr')
  maxadp = maxadd
  minadp = minadd
  alphao = 0e0
  alpha = 0e0
  beta = 0e0
!
! initialize the stopping test.
!
  call inithv (0)
  zthave = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
  itri = 1
  ip = itri
  if ( .not. (maxadd .or. minadd) ) go to 850
  ip = itri + 2*itmax
  call vfill (2*itmax,wk(itri),0e0)
 850  ir = ip + n
  iv1 = ir + n
  iv2 = iv1 + n
  nwusd = max(nwusd,iv2-1+n)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(ir))
!
!  begin iteration loop.
!
! determine whether or not to stop --
!
 10   call inithv (1)
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,wk(ir),wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
  if (in /= 0) go to 110
!
! perform first-iterate calculations
!
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(ir),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(ip))
  ard = vdot (n,wk(ip),wk(ip))
  go to 111
!
! perform subsequent-iterate calculations
!
 110  ardold = ard
!     if (abs(ardold) < srelpr) go to 996
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(ir),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  ard = vdot (n,wk(iv2),wk(iv2))
  an = ard/ardold
  call vtriad (n,wk(ip),wk(iv2),an,wk(ip),1)
  beta = an
!
! proceed to form the iterate.
!
 111  call suba (coef,jcoef,wfac,jwfac,n,wk(ip),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  pap = vdot (n,wk(iv2),wk(iv2))
  if (abs(pap) < srelpr**2) go to 998
  vlamda = ard/pap
!
  call vtriad (n,u,u,vlamda,wk(ip),1)
  call vtriad (n,wk(ir),wk(ir),-vlamda,wk(iv2),1)
!
! update eigenvalue estimates
!
  alphao = alpha
  alpha = vlamda
  if (maxadp .or. minadp) call chgcon (wk(itri),ier)
  if (ier < 0) go to 725
!
! proceed to next iteration
!
  in = in + 1
  is = is + 1
  go to 10
!
!  finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'cgnrw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' cgnr converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 995  ier = -16
  call ershow (ier,'cgnrw')
  return
!
 996  ier = -13
  call ershow (ier,'cgnrw')
  go to 725
!
 997  call ershow (ier,'cgnrw')
  go to 735
!
 998  ier = -15
  call ershow (ier,'cgnrw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'cgnrw')
  go to 735
!
end
subroutine cgw (suba,subq,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs,wksp,nw, &
  iparm,rparm,ier)
!
!*******************************************************************************
!
!! CGW drives the conjugate gradient algorithm.
!
!
!  Parameters:
!
!          suba   matrix-vector multiplication routine
!          subq   preconditioning routine
!          n      input integer.  order of the system (= nn)
!          u      input/output vector.  on input, u contains the
!                 initial guess to the solution.  on output, it
!                 contains the latest estimate to the solution.
!          ubar   input vector containing the true solution
!                  (optional)
!          rhs    input vector.  contains the right hand side
!                 of the matrix problem.
!          wksp   vector used for working space.
!          nw     length of wksp array.  if this length is less than
!                  the amount needed, nw will give the needed amount
!                  upon output.
!          iparm  integer vector of length 30.  allows user to
!                 specify some integer parameters which affect
!                 the method.
!          rparm  real vector of length 30.  allows user to
!                 specify some real parameters which affect
!                 the method.
!          ier    output integer.  error flag.
!
!  
!
  external  suba, subq
  integer   iparm(30), jcoef(2), jwfac(1)
  dimension rhs(1), u(1), ubar(1), wksp(1), rparm(30), coef(1), wfac(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, &
        iplr, iqlr, ntest, is, iacel, idgts, nbl1d, nbl2d
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, sigma, &
         rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
!  initialize common blocks
!
  ier = 0
  n = nn
  t1 = timer (dummy)
  iacel = 1
  timit = 0.0
  digit1 = 0.0
  digit2 = 0.0
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 35
  if (level >= 2) write (nout,10)
 10   format (1x,'cg')
!
!  compute workspace base addresses and check for sufficient
!  workspace.
!
  iw1 = 1
  iw2 = iw1 + n
  iw3 = iw2 + n
  iw4 = iw3 + n
  nwksp = 3*n + 2*itmax
  if (nw >= nwksp) go to 15
  ier = -2
  call ershow (ier,'cgw')
  go to 30
 15   continue
  call nmcalc (coef,jcoef,wfac,jwfac,1,subq,n,rhs,ubar,wksp,ier)
  if (ier < 0) go to 30
!
!  zero out workspace
!
  call vfill (nwksp,wksp,0.0)
!
!  iteration sequence
!
  call itcg (suba,subq,coef,jcoef,wfac,jwfac,n,u,ubar,rhs,wksp(iw1), &
    wksp(iw2),wksp(iw3),wksp(iw4),ier)
!
  if (ier < 0  .or.  ier == 1) go to 25
!
!  method has converged
!
  if (level >= 1) write (nout,20) in
 20   format (/1x,'cg  has converged in ',i5,' iterations' )
!
!  optional error analysis
!
 25   if (idgts < 0) go to 30
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wksp,digit1,digit2,idgts)
!
!  set return parameters in iparm and rparm
!
 30   t2 = timer (dummy)
  nw = 3*n + 2*in
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
  rparm(9) = omega
  rparm(10) = alphab
  rparm(11) = betab
  rparm(12) = specr
!
 35   continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
!
  return
end
subroutine chgcon (tri,ier)
!
!*******************************************************************************
!
!! CHGCON computes the new estimates for the largest and smallest eigenvalues.
!
!
!  Discussion:
!
!    These estimates are used for conjugate gradient acceleration.
!
!  Parameters:
!
!          tri    tridiagonal matrix associated with the eigenvalues
!                    of the conjugate gradient polynomial
!          ier    error code
!
!  
!
  dimension tri(2,2)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, &
         iplr, iqlr, ntest, is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, sigma, &
         rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
!     description of variables in common blocks in main routine
!
  save tl1,tl2,bl1,bl2
  ip = is
  if (ip - 1) 10,20,30
!
!  ip = 0
!
 10 end = 1.0/alpha
  tri(1,1) = end
  tri(2,1)= 0.0
  if (maxadp) emax = end
  if (minadp) emin = end
  return
!
!  ip = 1
!
 20   t1 = 1.0/alpha + beta/alphao
  t2 = beta/(alphao**2)
  tri(1,2) = t1
  tri(2,2) = t2
  tsqr = sqrt (t2)
  tl1 = tri(1,1) + tsqr
  tl2 = t1 + tsqr
  bl1 = tri(1,1) - tsqr
  bl2 = t1 - tsqr
  t3 = tri(1,1) + t1
  t4 = sqrt ( (t1-tri(1,1))**2 + 4.0*t2 )
  if (maxadp) emax = (t3 + t4)/2.0
  if (minadp) emin = (t3 - t4)/2.0
  return
!
!  ip >= 2
!
 30   t1 = 1.0/alpha + beta/alphao
  t2 = beta/(alphao**2)
  tsqr = sqrt (t2)
  tri(1,ip+1) = t1
  tri(2,ip+1) = t2
  if (.not. maxadp) go to 40
!
!  compute new estimate of emax.
!
  tl1 = amax1 (tl1,tl2+tsqr)
  tl2 = t1 + tsqr
  emaxo = emax
end = amax1 (tl1,tl2)
  e1 = eigvss (ip+1,tri,emaxo,end,2,ier)
  if (ier /= 3  .and.  ier /= 4) go to 35
!
!  poor estimate for emax.  therefore need to stop adaptive
!     procedure and keep old value of emax.
!
  maxadp = .false.
  if (level >= 2) write (nout,31) ier,in,emaxo
 31   format (/5x,'estimation of maximum eigenvalue emax halted' &
             /5x,'routine zbrent returned ier = ',i5 &
             /5x,'adaptive procedure turned off at iteration ',i5 &
             /5x,'final estimate of maximum eigenvalue =',e15.7/)
  go to 40
!
!  valid emax estimate.  check for small relative change in emax.
!
 35   emax = e1
  if (abs (emax - emaxo) < emax*zeta) maxadp = .false.
!
!  compute new estimate of emin.
!
 40   if (.not. minadp) return
  bl1 = amin1 (bl1,bl2-tsqr)
  bl2 = t1 - tsqr
  start = amax1 ( 0.0, amin1 (bl1,bl2) )
  emino = emin
  e1 = eigvss (ip+1,tri,start,emino,1,ier)
  if (ier /= 3  .and.  ier /= 4) go to 45
!
!  poor estimate for emin.  therefore need to stop adaptive
!     procedure and keep old value of emin.
!
  minadp = .false.
  if (level >= 2) write (nout,41) ier,in,emino
 41   format (/5x,'estimation of minimum eigenvalue emin halted' &
             /5x,'routine zbrent returned ier = ',i5 &
             /5x,'adaptive procedure turned off at iteration ',i5 &
             /5x,'final estimate of minimum eigenvalue =',e15.7/)
  return
!
!  valid emin estimate.  check for small relative change in emin.
!
 45   emin = e1
  if (abs (emin - emino) < emin*zeta) minadp = .false.
  return
end
subroutine chgsi (suba,coef,jcoef,wfac,jwfac,nn,z,wksp,icode,ier)
!
!*******************************************************************************
!
!! CHGSI adapts on the iteration parameters.
!
!
!  Parameters:
!
!         n         order of system (= nn)
!         z         current pseudo-residual vector
!         wksp      workspace vector of length n
!         icode     output indicator of parameter changes
!                    = 0    estimates of emax, emin not changed
!                    = 1    estimates of emax, emin changed
!         ier       error code
!
!  
!
  external suba
  integer jcoef(2), jwfac(1)
  dimension z(1), wksp(1), coef(1), wfac(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, &
        iplr, iqlr, ntest, is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
  n = nn
!
  istar = 3
  icode = 0
  if (is == 0) return
  rnrm = sqrt (rzdot)
  rnrmq = sqrt (dkq)
  rnrm1 = sqrt (dkm1)
  qa = rnrm/rnrmq
  t1 = rr**is
  qt = 2.0*sqrt (t1)/(1.0 + t1)
  if (qa <= qt**ff) return
  if (qa <= 1.0  .and.  is <= istar) return
  icode = 1
!
!  compute rayleigh quotient.
!         rq = (z,a*z)/(r,z)
!
  call suba (coef,jcoef,wfac,jwfac,n,z,wksp)
  top= 0.0
  do 10 i = 1,n
 10   top = top + z(i)*wksp(i)
  if (top >= 0.0) go to 15
  ier = -6
  call ershow (ier,'chgsi')
  return
 15   rq = top/rzdot
  kode = 0
  if (rq > rqmax) kode = 1
  rqmin = amin1 (rq,rqmin)
  rqmax = amax1 (rq,rqmax)
  yy = (1.0+t1)*(qa+sqrt (qa*qa-qt*qt))/2.0
  xx = yy**(1.0/float (is))
  if (qa > 1.0) go to 25
  if (kode == 1) go to 25
!
!  emin adjustment.
!
  eminp = (emax+emin)*(1.0-xx)*(xx-rr)/(2.0*xx*(rr+1.0))
  if (minadp) emin = amin1 (emin,eminp,rqmin)
  if (maxadp) emax = amax1 (emax,rqmax)
  if (level >= 2) write (nout,20) in,rq,eminp,emin,emax
 20   format (/1x,15x,'parameters were changed at iteration',i7/ &
      1x,20x,'rayleigh quotient  ',f15.9/ &
      1x,20x,'young estimate     ',f15.9/ &
      1x,20x,'emin               ',f15.9/ &
      1x,20x,'emax               ',f15.9/)
  return
!
!  emax adjustment.
!
 25   emaxp = (emax+emin)*(1.0+xx)*(xx+rr)/(2.0*xx*(rr+1.0))
  uu = ((1.0+t1)/(1.0+rr**(is-1))) * (rnrm/rnrm1)
  emaxpp = (emax+emin)*(1.0+uu)*(uu+rr)/(2.0*uu*(rr+1.0))
  if (maxadp) emax = amax1 (emax,1.1*emaxp,1.1*emaxpp,1.1*rqmax)
  if (minadp) emin = rqmin
  if (level >= 2) write (nout,30) in,rq,emaxp,emaxpp,emin,emax
 30   format (/1x,15x,'parameters were changed at iteration',i7/ &
     20x,'rayleigh quotient  ',f15.9/ &
     20x,'young estimate     ',f15.9/ &
     20x,'hageman estimate   ',f15.9/ &
     20x,'emin               ',f15.9/ &
     20x,'emax               ',f15.9/)
  return
end
subroutine ckconv (ier)
!
!*******************************************************************************
!
!! CKCONV checks if the iterative method has stagnated or had other misfortunes.
!
  parameter (nst=20)
  parameter (eps=1.e-7)
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
!
  dimension stold(nst)
  save stold, ist
  ind(i) = 1 + mod(i,nst)
!
  if (in <= 0) ist = 0
!
  ist = ist + 1
  stold(ind(ist)) = stptst
  if (ist < nst) go to 900
!     do 2 i = 1, nst-1
  do 2 i = nst-1, 1, -1
!     val = abs(stold(ind(ist-i))-stptst)/stptst
  val = abs(stold(ind(ist-i))-stptst)
  if (val > eps*stptst) go to 900
 2    continue
  ier = -19
  call ershow (ier,'ckconv')
  return

 900  return
end
subroutine color ( nxp, nyp, nzp, nx, ny, nz, pp, p )
!
!*******************************************************************************
!
!! COLOR expands a color pattern to a full grid color array.
!
!
!  Discussion:
!
!    The (small) color pattern array PP is repeatedly mapped onto
!    the large grid color array P, in the same way that a 2 by 2
!    block of red/black squares can be used to define the color
!    pattern on a large checkerboard.
!
!  Parameters:
!
!       nxp,    integer variables giving the x, y, and z dimensions
!        nyp,    of the pattern array, respectively.
!        nzp
!       nx,ny,  integer variables giving the x, y, and z dimensions
!        nz      of the grid, respectively.
!       pp      integer vector of length  nxp*nyp*nzp
!                giving the color pattern to be repeated
!       p       integer vector of length  nxg*nyg*nzg
!                which contains upon output the grid coloring
!
  integer nx
  integer nxp
  integer ny
  integer nyp
  integer nz
  integer nzp
!
  integer i
  integer ip
  integer j
  integer jp
  integer k
  integer kp
  integer p(nx,ny,nz)
  integer pp(nxp,nyp,nzp)
!
  do k = 1, nz

    kp = mod ( k - 1, nzp ) + 1

    do j = 1, ny

      jp = mod ( j - 1, nyp ) + 1

      do i = 1, nx

        ip = mod ( i - 1, nxp ) + 1

        p(i,j,k) = pp(ip,jp,kp)

      end do
    end do
  end do

  return
end
subroutine copy (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! COPY does a vector copy (null preconditioner).
!
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  do 10 i = 1,n
 10   z(i) = r(i)
  return
end
subroutine defcon (ndim,nn,maxnz,jcoef,coef,kblsz,iblock,lbhb)
!
!*******************************************************************************
!
!! DEFCON defines block constants for block-structured matrices.
!
!
!     (diagonal data structure, constant block size)
!
!  Parameters:
!
!         ndim     row dimension of coef array in defining routine
!         nn       size of system
!         maxnz    number of diagonals in coef
!         jcoef    integer vector of size maxnz giving the diagonal
!                   numbers
!         coef     matrix representation array
!         kblsz    constant block size
!         iblock   integer array of size 3 by lbhb
!                   giving block constants upon output
!         lbhb     integer giving the number of diagonal blocks
!                   upon output.
!
!  
!
  integer   jcoef(2), iblock(3,3)
  dimension coef(ndim,1)
!
  n = nn
  ipt = 2
  iblock(1,1) = 0
  iblock(1,2) = 0
  iblock(2,1) = 1
  iblock(3,1) = 0
  iblock(3,2) = 0
  do 25 j = 1,maxnz
     jd = jcoef(j)
     do 10 i = 1,n
        if (coef(i,j) /= 0.0) go to 15
 10      continue
     go to 25
 15      jcol = i + jd
!
!  find block for jcol.
!
     ib = (i-1)/kblsz + 1
     jb = (jcol-1)/kblsz + 1
     id = jb - ib
     if (id == iblock(1,ipt)) go to 20
     ipt = ipt + 1
     iblock(1,ipt) = id
     iblock(3,ipt) = 0
 20      iblock(3,ipt) = iblock(3,ipt) + 1
 25   continue
  lbhb = ipt
!
!  split zero diagonal block into super and sub diagonals.
!
  jlim = iblock(3,2)
  do 30 j = 1,jlim
     jd = jcoef(j)
     if (jd < 0) go to 35
     iblock(3,1) = iblock(3,1) + 1
     iblock(3,2) = iblock(3,2) - 1
 30   continue
  j = jlim + 1
 35   iblock(2,2) = j
!
!  form starting positions.
!
  if (lbhb <= 2) return
  iblock(2,3) = 1
  if (lbhb <= 3) return
  do 40 j = 4,lbhb
 40   iblock(2,j) = iblock(2,j-1) + iblock(3,j-1)
  return
end
subroutine define (ndim,maxnew,jcnew,coef,ncol,nc,iblock,lbhb)
!
!*******************************************************************************
!
!! DEFINE defines block constants for block-structured matrices.
!
!
!     (diagonal data structure, nonconstant block size)
!
!  Parameters:
!
!         ndim     row dimension of coef array in defining routine
!         maxnew   integer vector giving the number of diagonals
!                   for each distinct block size.
!         jcnew    integer array of size ncolor*max(maxnew(i))
!                   giving the diagonal numbers for each distinct
!                   block size.
!         coef     matrix representation array
!         ncolor   number of distinct block sizes
!         nc       integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants upon output
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size
!                   upon output.
!
!  
!
  integer   maxnew(ncol), jcnew(ncol,1), nc(ncol), lbhb(ncol), iblock(3,ncol,3)
  dimension coef(ndim,1)
!
  ncolor = ncol
  ist = 1
  do 60 k = 1,ncolor
     ncc = nc(k)
     maxnz = maxnew(k)
     ied = ist + ncc - 1
     ipt = 2
     iblock(1,k,1) = 0
     iblock(1,k,2) = 0
     iblock(2,k,1) = 1
     iblock(3,k,1) = 0
     iblock(3,k,2) = 0
     do 35 j = 1,maxnz
        jd = jcnew(k,j)
        do 10 i = ist,ied
           if (coef(i,j) /= 0.0) go to 15
 10         continue
        go to 35
 15         jcol = i + jd
!
!  find block for jcol.
!
        ib = k
        js = 0
        do 20 ij = 1,ncolor
           js = js + nc(ij)
           if (js >= jcol) go to 25
 20         continue
 25         jb = ij
        id = jb - ib
        if (id == iblock(1,k,ipt)) go to 30
        ipt = ipt + 1
        iblock(1,k,ipt) = id
        iblock(3,k,ipt) = 0
 30         iblock(3,k,ipt) = iblock(3,k,ipt) + 1
 35      continue
     lbhb(k) = ipt
!
!  split zero diagonal block into super and sub diagonals.
!
     jlim = iblock(3,k,2)
     do 40 j = 1,jlim
        jd = jcnew(k,j)
        if (jd < 0) go to 45
        iblock(3,k,1) = iblock(3,k,1) + 1
        iblock(3,k,2) = iblock(3,k,2) - 1
 40      continue
     j = jlim + 1
 45      iblock(2,k,2) = j
!
!  form starting positions.
!
     jlim = lbhb(k)
     if (jlim <= 2) go to 55
     iblock(2,k,3) = 1
     if (jlim <= 3) go to 55
     do 50 j = 4,jlim
 50      iblock(2,k,j) = iblock(2,k,j-1) + iblock(3,k,j-1)
 55      ist = ied + 1
 60   continue
  return
end
function determ (n,tri,xlmda)
!
!*******************************************************************************
!
!! DETERM computes the determinant of a symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    The matrix is given by tri. det(tri - xlmda*i) = 0
!
!  Parameters:
!
!          n      order of tridiagonal system
!          tri    symmetric tridiagonal matrix of order n
!          xlmda  argument for characteristic equation
!
!  
!
  dimension tri(2,1)
!
  nm1 = n - 1
  d2 = tri(1,n) - xlmda
  d1 = d2 * (tri(1,nm1) - xlmda) - tri(2,n)
     if (n == 2) go to 20
!
!  beginning of loop
!
  do 10 l = nm1,2,-1
     d3 = d2
     d2 = d1
     d1 = (tri(1,l-1) - xlmda) * d2 - d3 * tri(2,l)
   10 continue
!
!  determinant computed
!
   20 determ = d1
!
  return
end
subroutine detsym (ndim,maxnzz,coef,jcoef,nn,isymm)
!
!*******************************************************************************
!
!! DETSYM determines if the matrix is symmetric.
!
!
!     (Purdue storage format)
!
!  Parameters:
!
!         ndim     row dimension of coef in defining routine
!         maxnz    number of columns in coef
!         coef     array of matrix nonzeros
!         jcoef    array of matrix column numbers
!         n        order of matrix (= nn)
!         isymm    symmetry switch.  upon output,
!                   isymm = 0  if matrix is symmetric
!                         = 1  if matrix is nonsymmetric
!
!  
!
  dimension coef(ndim,2)
  integer   jcoef(ndim,2)
!
  n = nn
  maxnz = maxnzz
  isymm = 0
  if (maxnz <= 1) return
  do 20 i = 1,n
     do 15 j = 2,maxnz
        jcol = jcoef(i,j)
        if (jcol == i) go to 15
        val = coef(i,j)
        do 10 jj = 2,maxnz
           jcol1 = jcoef(jcol,jj)
           if (jcol1 /= i) go to 10
           val1 = coef(jcol,jj)
           if (val1 == val) go to 15
           isymm = 1
           return
 10         continue
        isymm = 1
        return
 15      continue
 20   continue
  return
end
subroutine dfault ( iparm, rparm )
!
!*******************************************************************************
!
!! DFAULT sets the default values of IPARM and RPARM.
!
!
!  Parameters:
!
!          iparm
!           and
!          rparm  arrays specifying options and tolerances
!
!
!  
!
  integer   iparm(30)
  dimension rparm(30)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
!     description of variables in common blocks in main routine
!
!     srelpr  - computer precision (approx.)
!     if installer of package does not know srelpr value,
!     an approximate value can be determined from a simple
!     fortran program such as
!
!     srelpr = 1.0
!  2  srelpr = 0.5*srelpr
!     temp = srelpr + 1.0
!     if (temp > 1.0) go to 2
!     srelpr = 2.0*srelpr
!     write (6,3) srelpr
!  3  format (1x,'srelpr = ',e20.10)
!     stop
!   end
!
!
!     some values are-
!
!     srelpr = 7.1e-15   for cray x-mp  (approx.) 2**-47
!            = 1.49e-8   for dec 10  (approx.) 2**-26
!            = 1.192e-7  for vax 11/780 (approx) 2**-23
!            = 4.768e-7  for ibm 370/158
!
!             *** should be changed for other machines ***
!
!     to facilitate convergence, rparm(1) should be set to
!          500.*srelpr or larger
!
!
  srelpr = epsilon ( srelpr )
!
!  keygs is a flag to specify how gather/scatter operations
!     are performed.
!       = 1    gather explicitly into a workspace vector
!       = 2    gather implicitly using indirect addressing
!
!
  keygs = 1
!
!  keyzer is a flag to specify if memory has been zeroed out.
!     i.e., is the operation  0.0 * indefinite = 0.0  legal
!       = 0    not legal
!       = 1    legal
!
  keyzer = 0
!
  iparm(1)  =      2
  iparm(2)  =    100
  iparm(3)  =      0
  iparm(4)  =      6
  iparm(5)  =      0
  iparm(6)  =      1
  iparm(7)  =      1
  iparm(8)  =      1
  iparm(9)  =      5
  iparm(10) = 100000
  iparm(11) =      0
  iparm(12) =      2
  iparm(13) =      0
  iparm(14) =      0
  iparm(15) =      1
  iparm(16) =      0
  iparm(17) =      0
  iparm(18) =      2
  iparm(19) =     -1
  iparm(20) =     -1
  iparm(21) =      1
  iparm(22) =      1
  iparm(23) =      2
  iparm(24) =      0
  iparm(25) =      1
!
  rparm(1)  = 1.0e-4
  rparm(2)  = 2.0
  rparm(3)  = 1.0
  rparm(4)  = 0.75
  rparm(5)  = 0.75
  rparm(6)  = 0.0
  rparm(7)  = 0.0
  rparm(8)  = 0.0
  rparm(9)  = 1.0
  rparm(10) = 0.0
  rparm(11) = 0.25
  rparm(12) = 0.0
  rparm(13) = 0.0
  rparm(14) = 0.0
  rparm(15) = 500.0 * srelpr
  rparm(16) = 0.0
!
  return
end
subroutine echall (n,iparm,rparm,icall,icallr,ier)
!
!*******************************************************************************
!
!! ECHALL initializes the package common blocks.
!
!
!  It uses the information contained in iparm and rparm.  echall also
!  prints the values of all parameters in iparm and rparm.
!
!  Parameters:
!
!          iparm
!           and
!          rparm  arrays of parameters specifying options and
!                    tolerances
!          icall  indicator of which parameters are being printed
!                    icall = 1,  initial parameters
!                          = 2,  final parameters
!          icallr  indicator of calling routine
!                          = 1,  called from nspcg
!                          = 2,  called from accelerator
!
!  
!
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / itcom8 / ainf
!
!
!
  logical erflag
  integer   iparm(25)
  dimension rparm(16)
  character*6 inames(25), rnames(16)
  data naiprm, narprm / 11,12 /
  data inames / 'ntest', 'itmax', 'level', 'nout', 'idgts', &
    'maxadp', 'minadp', 'iomgad', 'ns1', 'ns2', 'ns3', &
    'nstore', 'iscale', 'iperm', 'ifact', 'lvfill', &
    'ltrunc', 'ipropa', 'kblsz', 'nbl2d', 'ifctv', &
    'iqlr', 'isymm', 'ielim', 'ndeg' /

  data rnames / 'zeta', 'emax', 'emin', 'ff', 'fff', 'timit', &
    'digit1', 'digit2', 'omega', 'alphab', 'betab', &
    'specr', 'timfac', 'timtot', 'tol', 'ainf' /
!
  if (icall /= 1) go to 20
!
! handle accelerator parameters.
!
  ntest  = iparm(1)
  itmax  = iparm(2)
  level  = iparm(3)
  nout   = iparm(4)
  idgts  = iparm(5)
  maxad  = iparm(6)
  maxadd = (maxad == 1)
  minad  = iparm(7)
  minadd = (minad == 1)
  maxadp = maxadd
  minadp = minadd
  iomgad = iparm(8)
  omgadp = (iomgad == 1)
  ns1    = iparm(9)
  ns2    = iparm(10)
  ns3    = iparm(11)
  iqlr   = iparm(22)
  iplr   = iqlr
!
  zeta   = rparm(1)
  emax   = rparm(2)
  emin   = rparm(3)
  ff     = rparm(4)
  fff    = rparm(5)
  timit  = rparm(6)
  digit1 = rparm(7)
  digit2 = rparm(8)
  omega  = rparm(9)
  alphab = rparm(10)
  betab  = rparm(11)
  specr  = rparm(12)
!
  erflag = .false.
  erflag = erflag .or. ntest < 1 .or. ntest > 10
  erflag = erflag .or. itmax <= 0
  erflag = erflag .or. maxad < 0 .or. maxad > 1
  erflag = erflag .or. minad < 0 .or. minad > 1
  erflag = erflag .or. ns1 < 0
  erflag = erflag .or. ns2 < 0
  erflag = erflag .or. emax < 0.0
  erflag = erflag .or. emin < 0.0
  erflag = erflag .or. ff <= 0.0 .or. ff > 1.0
  if (erflag) go to 999
!
!  test if eps is too small
!
  temp = 500.0*srelpr
  if (zeta >= temp) go to 150
  ier = 2
  call ershow (ier,'echall')
  zeta = temp
  rparm(1) = temp
!
!  verify n
!
 150  if (n > 0 ) go to 200
  ier = -1
  call ershow (ier,'echall')
  return
!
! now handle preconditioner parameters.
!
 200  if (icallr == 2) go to 50
  nstore = iparm(12)
  iscale = iparm(13)
  iperm  = iparm(14)
  ifact  = iparm(15)
  lvfill = iparm(16)
  ltrunc = iparm(17)
  ipropa = iparm(18)
  nbl1d  = iparm(19)
  nbl2d  = iparm(20)
  ifctv  = iparm(21)
  iqlr   = iparm(22)
  isymm  = iparm(23)
  ndeg   = iparm(25)
  ainf   = rparm(16)
!
  if (nbl1d == -1) nbl1d = n
  if (nbl2d == -1) nbl2d = n
  kblsz = nbl1d
  erflag = .false.
  erflag = erflag .or. iqlr < 0 .or. iqlr > 3
  erflag = erflag .or. ipropa < 0 .or. ipropa > 3
  if (erflag) go to 999
!
!
!  initialize rest of common variables
!
 50   halt   = .false.
  stptst= 0.0
  udnm   = 1.0
  in     = 0
!
!  Prepare to do output.
!
  if (level <= 2) return

  write ( nout, * ) ' '
  write ( nout, * ) 'Initial iterative parameters'
  write ( nout, * ) ' '

  go to 30

 20   if (level <= 2) return

  write ( nout, * ) ' '
  write ( nout, * ) 'Final iterative parameters'
  write ( nout, * ) ' '

 30   if (icallr == 2) go to 305

  write ( nout, * ) ' '
  write ( nout, * ) 'Preprocessor and preconditioner parameters'
  write ( nout, * ) ' '

  ibip = naiprm + 1
  ieip = 25
  ibrp = narprm + 1
  ierp = 16
  go to 300
 305  continue
  write ( nout, * ) ' '
  write ( nout, * ) 'General and acceleration parameters'
  write ( nout, * ) ' '

  ibip = 1
  ieip = naiprm
  ibrp = 1
  ierp = narprm

 300  write (nout,35) (i,iparm(i),inames(i),i=ibip,ieip)
 35   format (10x,'iparm(',i2,') =',i15,4x,'(',a6,')'  )
  write (nout,40) (i,rparm(i),rnames(i),i=ibrp,ierp)
 40   format (10x,'rparm(',i2,') =',e15.8,4x,'(',a6,')'  )
  return
!
! error returns.
!
! inadmissible option.
 999  ier = -10
  call ershow (ier,'echall')
  return
end
function eigvss (n,tri,start,end,icode,ier)
!
!*******************************************************************************
!
!! EIGVSS computes a selected eigenvalue of a symmetric tridiagonal matrix.
!
!
!  Discussion:
!
!    The eigenvalue is computed for conjugate gradient acceleration.
!    The modified imsl routine zbrent is used.
!
!  Parameters:
!
!          n      order of tridiagonal system
!          tri    symmetric tridiagonal matrix of order n
!          start  initial lower bound of interval containing root
!        end    initial upper bound of interval containing root
!          icode  operation key
!                   = 1   minimum eigenvalue sought
!                   = 2   maximum eigenvalue sought
!          ier    error flag
!
!  
!
  dimension tri(2,1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  eigvss= 0.0
  itmp = int (-alog10 (abs (zeta)))
  nsig = max (itmp,4)
  maxfn = max (itmax,50)
  eps= 0.0
  a = start
  b = end
  call zbrent (n,tri,eps,nsig,a,b,maxfn,ier)
  if (icode == 1) eigvss = amax1 (a,b)
  if (icode == 2) eigvss = amin1 (a,b)
!
  return
end
subroutine elim (n,jcoef,coef,rhs,wksp,iwksp,toll)
!
!*******************************************************************************
!
!! ELIM removes certains rows of the matrix.
!
!
!  Discussion:
!
!    The eliminated rows are those for which the ratio of the
!    sum of off-diagonal elements to the diagonal element is
!    small (less than tol) in absolute value.
!
!    this is to take care of matrices arising from finite
!    element discretizations of partial differential equations
!    with dirichlet boundary conditions implemented by penalty
!    methods.  any such rows and corresponding columns are then
!    eliminated (set to the identity after correcting the rhs).
!
!  parameter list --
!
!         n       dimension of matrix
!         jcoef   integer array of matrix representation
!         coef    array of sparse matrix representation
!         rhs     right hand side of matrix problem
!         wksp    wksp array of length n
!         tol     tolerance factor  (= toll)
!
!  specifications for arguments
!
  common / cmpart / mpstrt, mpart
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer   jcoef(2), iwksp(1)
  dimension coef(1), rhs(1), wksp(1)
!
  tol = toll
  go to (5,10,15,20,25), nstore
 5    call elim1 (n,ndim,maxnz,jcoef,coef,rhs,wksp(irpnt),tol)
  return
 10   call elim2 (n,ndim,maxnz,jcoef,coef,rhs,wksp(irpnt),tol)
  return
 15   call elim3 (n,ndim,maxnz,jcoef,coef,rhs,wksp(irpnt),tol)
  return
 20   continue

  call elim4 (mpart,iwksp(mpstrt),jcoef,jcoef(ndim+1),coef,rhs,wksp(irpnt),tol)
  return
 25   continue

  call elim5 (mpart,iwksp(mpstrt),jcoef,jcoef(ndim+1),coef,rhs,wksp(irpnt),tol)
  return
end
subroutine elim1 (nn,ndim,maxnzz,jcoef,coef,rhs,wksp,toll)
!
!*******************************************************************************
!
!! ELIM1 removes certina rows of the matrix.
!
!
!  Discussion:
!
!    The elminated rows are those for which the ratio of the
!    sum of off-diagonal elements to the diagonal element is
!    small (less than tol) in absolute value.
!
!    this is to take care of matrices arising from finite
!    element discretizations of partial differential equations
!    with dirichlet boundary conditions implemented by penalty
!    methods.  any such rows and corresponding columns are then
!    eliminated (set to the identity after correcting the rhs).
!    Purdue format.
!
!  parameter list --
!
!         n       dimension of matrix ( = nn)
!         ndim    row dimension of arrays jcoef and coef in the
!                    calling program
!         maxnz   maximum number of nonzero entries per row (=maxnzz)
!         jcoef   integer array of matrix representation
!         coef    array of sparse matrix representation
!         rhs     right hand side of matrix problem
!         wksp    wksp array of length n
!         tol     tolerance factor  (= toll)
!
!  specifications for arguments
!
  integer   jcoef(ndim,1)
  dimension coef(ndim,1), rhs(1), wksp(1)
!
  n = nn
  maxnz = maxnzz
  tol = toll
  if (n <= 0  .or.  maxnz < 2) return
!
!  find maximum off-diagonal elements in absolute value.
!
  do 10 i = 1,n
 10   wksp(i)= 0.0
  do 20 j = 2,maxnz
     do 15 i = 1,n
 15      wksp(i) = wksp(i) + abs (coef(i,j))
 20   continue
  do 25 i = 1,n
 25   wksp(i) = wksp(i) / abs(coef(i,1))
!
!  eliminate desired rows and columns.
!
  do 35 i = 1,n
     if (wksp(i) > tol) go to 35
     rhs(i) = rhs(i)/coef(i,1)
     coef(i,1) = 1.0
     do 30 j = 2,maxnz
        coef(i,j)= 0.0
        jcoef(i,j) = i
 30      continue
 35   continue
  do 45 j = 2,maxnz
     do 40 i = 1,n
        jcol = jcoef(i,j)
        if (wksp(jcol) > tol) go to 40
        rhs(i) = rhs(i) - coef(i,j)*rhs(jcol)
        coef(i,j)= 0.0
        jcoef(i,j) = i
 40      continue
 45   continue
  return
end
subroutine elim2 (nn,ndim,maxnzz,jcoef,coef,rhs,wksp,toll)
!
!*******************************************************************************
!
!! ELIM2 removes certain rows of the matrix.
!
!
!  The eliminated rows are those for which the ratio of the
!     sum of off-diagonal elements to the diagonal element is
!     small (less than tol) in absolute value.
!     this is to take care of matrices arising from finite
!     element discretizations of partial differential equations
!     with dirichlet boundary conditions implemented by penalty
!     methods.  any such rows and corresponding columns are then
!     eliminated (set to the identity after correcting the rhs).
!     symmetric diagonal format.
!
!  parameter list --
!
!         n       dimension of matrix ( = nn)
!         ndim    row dimension of array coef in the
!                    calling program
!         maxnz   number of diagonals stored
!         jcoef   integer vector of diagonal numbers
!         coef    array of sparse matrix representation
!         rhs     right hand side of matrix problem
!         wksp    wksp array of length n
!         tol     tolerance factor  (= toll)
!
!  specifications for arguments
!
  integer   jcoef(1)
  dimension coef(ndim,1), rhs(1), wksp(1)
!
  n = nn
  maxnz = maxnzz
  tol = toll
  if (n <= 0  .or.  maxnz < 2) return
!
!  find maximum off-diagonal elements in absolute value.
!
  do 10 i = 1,n
 10   wksp(i)= 0.0
  do 25 j = 2,maxnz
     ind = jcoef(j)
     len = n - ind
     do 15 i = 1,len
 15      wksp(i) = wksp(i) + abs (coef(i,j))
     do 20 i = 1,len
 20      wksp(i+ind) = wksp(i+ind) + abs (coef(i,j))
 25   continue
  do 30 i = 1,n
 30   wksp(i) = wksp(i) / abs(coef(i,1))
!
!  eliminate desired rows and columns.
!
  do 50 i = 1,n
     if (wksp(i) > tol) go to 50
     rhs(i) = rhs(i)/coef(i,1)
     coef(i,1) = 1.0
     do 40 j = 2,maxnz
        jcol = jcoef(j)
        iback = i - jcol
        iforw = i + jcol
        if (iforw > n) go to 35
        if (wksp(iforw) <= tol) go to 35
        rhs(iforw) = rhs(iforw) - coef(i,j)*rhs(i)
 35         if (iback < 1) go to 40
        rhs(iback) = rhs(iback) - coef(iback,j)*rhs(i)
        coef(iback,j)= 0.0
 40      continue
     do 45 j = 2,maxnz
 45      coef(i,j)= 0.0
 50   continue
  return
end
subroutine elim3 (nn,ndim,maxnzz,jcoef,coef,rhs,wksp,toll)
!
!*******************************************************************************
!
!! ELIM3 removes certain rows of the matrix.
!
!
!  The eliminated rows are those for which the ratio of the
!     sum of off-diagonal elements to the diagonal element is
!     small (less than tol) in absolute value.
!     this is to take care of matrices arising from finite
!     element discretizations of partial differential equations
!     with dirichlet boundary conditions implemented by penalty
!     methods.  any such rows and corresponding columns are then
!     eliminated (set to the identity after correcting the rhs).
!     nonsymmetric diagonal format.
!
!  parameter list --
!
!         n       dimension of matrix ( = nn)
!         ndim    row dimension of array coef in the
!                    calling program
!         maxnz   number of diagonals stored
!         jcoef   integer vector of diagonal numbers
!         coef    array of sparse matrix representation
!         rhs     right hand side of matrix problem
!         wksp    wksp array of length n
!         tol     tolerance factor  (= toll)
!
!  specifications for arguments
!
  integer   jcoef(1)
  dimension coef(ndim,1), rhs(1), wksp(1)
!
  n = nn
  maxnz = maxnzz
  tol = toll
  if (n <= 0  .or.  maxnz < 2) return
!
!  find maximum off-diagonal elements in absolute value.
!
  do 10 i = 1,n
 10   wksp(i)= 0.0
  do 20 j = 2,maxnz
     ind = jcoef(j)
     ist1 = max(1,1 - ind)
     ist2 = min(n,n - ind)
     do 15 i = ist1,ist2
 15      wksp(i) = wksp(i) + abs (coef(i,j))
 20   continue
  do 25 i = 1,n
 25   wksp(i) = wksp(i) / abs(coef(i,1))
!
!  eliminate desired rows and columns.
!
  do 35 i = 1,n
     if (wksp(i) > tol) go to 35
     rhs(i) = rhs(i)/coef(i,1)
     coef(i,1) = 1.0
     do 30 j = 2,maxnz
 30      coef(i,j)= 0.0
 35   continue
  do 45 i = 1,n
     if (wksp(i) > tol) go to 45
     do 40 j = 2,maxnz
        inew = i - jcoef(j)
        if (inew < 1 .or. inew > n) go to 40
        rhs(inew) = rhs(inew) - coef(inew,j)*rhs(i)
        coef(inew,j)= 0.0
 40      continue
 45   continue
  return
end
subroutine elim4 (mm,np,ia,ja,a,rhs,wksp,toll)
!
!*******************************************************************************
!
!! ELIM4 removes certain rows of the matrix.
!
!
!  The eliminated rows are those for which the ratio of the
!     sum of off-diagonal elements to the diagonal element is
!     small (less than tol) in absolute value.
!     this is to take care of matrices arising from finite
!     element discretizations of partial differential equations
!     with dirichlet boundary conditions implemented by penalty
!     methods.  any such rows and corresponding columns are then
!     eliminated (set to the identity after correcting the rhs).
!     symmetric sparse format.
!
!  parameter list --
!
!         m       number of partitions
!         np      pointer vector to partitions
!         ia      vector of i values
!         ja      vector of j values
!         a       vector of coefficients
!         rhs     right hand side of matrix problem
!         wksp    wksp vector of length n (2n if keygs = 1)
!         tol     tolerance factor  (= toll)
!
!  specifications for arguments
!
  integer ia(1), ja(1), np(2)
  dimension a(1), rhs(1), wksp(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  m = mm
  n = np(2) - 1
  nz = np(m+1) - 1
  tol = toll
  np1 = n + 1
!
!  find sum of absolute values of off-diagonal coefficients.
!
  do 10 i = 1,n
 10   wksp(i)= 0.0
  if (keygs == 1) go to 30
  do 25 k = 2,m
     ist = np(k)
     ied = np(k+1) - 1
!dir$ ivdep
     do 15 i = ist,ied
 15      wksp(ia(i)) = wksp(ia(i)) + abs(a(i))
!dir$ ivdep
     do 20 i = ist,ied
 20      wksp(ja(i)) = wksp(ja(i)) + abs(a(i))
 25   continue
  go to 50
 30   do 45 k = 2,m
     ist = np(k)
     ied = np(k+1) - 1
     len = ied - ist + 1
     call vgathr (len,wksp,ia(ist),wksp(n+1))
     do 35 i = ist,ied
 35      wksp(i-ist+1+n) = wksp(i-ist+1+n) + abs(a(i))
     call vscatr (len,wksp(n+1),ia(ist),wksp)
     call vgathr (len,wksp,ja(ist),wksp(n+1))
     do 40 i = ist,ied
 40      wksp(i-ist+1+n) = wksp(i-ist+1+n) + abs(a(i))
     call vscatr (len,wksp(n+1),ja(ist),wksp)
 45   continue
 50   do 55 i = 1,n
 55   wksp(i) = wksp(i) / abs(a(i))
!
!  eliminate desired rows and columns.
!
  do 70 l = 1,n
     if (wksp(l) > tol) go to 70
     rhs(l) = rhs(l)/a(l)
     a(l) = 1.0
     do 60 k = np1,nz
        i = ia(k)
        j = ja(k)
        if (i == l .and. wksp(j) > tol) rhs(j) = rhs(j) - a(k)*rhs(i)
        if (j /= l) go to 60
        rhs(i) = rhs(i) - a(k)*rhs(j)
        a(k) = 0.0
 60      continue
     do 65 k = np1,nz
        if (ia(k) == l) a(k) = 0.0
 65      continue
 70   continue
  return
end
subroutine elim5 (mm,np,ia,ja,a,rhs,wksp,toll)
!
!*******************************************************************************
!
!! ELIM5 removes certain rows of the matrix.
!
!  The elminated rows are those for which the ratio of the
!     sum of off-diagonal elements to the diagonal element is
!     small (less than tol) in absolute value.
!     this is to take care of matrices arising from finite
!     element discretizations of partial differential equations
!     with dirichlet boundary conditions implemented by penalty
!     methods.  any such rows and corresponding columns are then
!     eliminated (set to the identity after correcting the rhs).
!     nonsymmetric sparse format.
!
!  parameter list --
!
!         m       number of partitions
!         np      pointer vector to partitions
!         ia      vector of i values
!         ja      vector of j values
!         a       vector of coefficients
!         rhs     right hand side of matrix problem
!         wksp    wksp vector of length n (2n if keygs = 1)
!         tol     tolerance factor  (= toll)
!
!  specifications for arguments
!
  integer ia(1), ja(1), np(2)
  dimension a(1), rhs(1), wksp(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  m = mm
  n = np(2) - 1
  nz = np(m+1) - 1
  tol = toll
!
!  find sum of absolute values of off-diagonal coefficients.
!
  do 10 i = 1,n
 10   wksp(i)= 0.0
  if (keygs == 1) go to 25
  do 20 k = 2,m
     ist = np(k)
     ied = np(k+1) - 1
!dir$ ivdep
     do 15 i = ist,ied
 15      wksp(ia(i)) = wksp(ia(i)) + abs(a(i))
 20   continue
  go to 40
 25   do 35 k = 2,m
     ist = np(k)
     ied = np(k+1) - 1
     len = ied - ist + 1
     call vgathr (len,wksp,ia(ist),wksp(n+1))
     do 30 i = ist,ied
 30      wksp(i-ist+1+n) = wksp(i-ist+1+n) + abs(a(i))
     call vscatr (len,wksp(n+1),ia(ist),wksp)
 35   continue
 40   do 45 i = 1,n
 45   wksp(i) = wksp(i) / abs(a(i))
!
!  eliminate desired rows and columns.
!
  do 50 i = 1,n
     if (wksp(i) > tol) go to 50
     rhs(i) = rhs(i)/a(i)
     a(i) = 1.0
 50   continue
  np1 = n + 1
  do 55 k = np1,nz
     if (wksp(ia(k)) <= tol) a(k) = 0.0
 55   continue
  do 60 k = np1,nz
     j = ja(k)
     if (wksp(j) > tol) go to 60
     i = ia(k)
     rhs(i) = rhs(i) - a(k)*rhs(j)
     a(k) = 0.0
 60   continue
  return
end
subroutine ershow (ierr,iname)
!
!*******************************************************************************
!
!! ERSHOW prints an appropriate error message for the error numbered IER.
!
!
!  Parameters:
!
!        ier     error number (input)
!                 > 0     for warning errors
!                 < 0     for fatal errors
!        iname   routine name in which error occurred
!
!  
!
  character*(*) iname
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  character*80  fmess(20), wmess(6)
  data fmess(1)  / 'nonpositive matrix size n' /
  data fmess(2)  / 'insufficient real workspace' /
  data fmess(3)  / 'insufficient integer workspace' /
  data fmess(4)  / 'nonpositive diagonal element' /
  data fmess(5)  / 'nonexistent diagonal element' /
  data fmess(6)  / 'a is not positive definite' /
  data fmess(7)  / 'q is not positive definite' /
  data fmess(8)  / 'unable to permute matrix as requested' /
  data fmess(9)  / 'mdim not large enough to expand matrix' /
  data fmess(10) / 'inadmissible parameter encountered' /
  data fmess(11) / 'incorrect storage mode for block method' /
  data fmess(12) / 'zero pivot encountered during factorization' /
  data fmess(13) / 'breakdown in direction vector calculation' /
  data fmess(14) / 'breakdown in attempt to perform rotation' /
  data fmess(15) / 'breakdown in iterate calculation' /
  data fmess(16) / 'unimplemented combination of parameters' /
  data fmess(17) / 'error in computing preconditioning polynomial' /
  data fmess(18) / 'unable to perform eigenvalue estimation' /
  data fmess(19) / 'iterative method has gone to sleep' /
  data fmess(20) / 'unknown error' /
  data wmess(1)  / 'failure to converge in itmax iterations' /
  data wmess(2)  / 'zeta too small' /
  data wmess(3)  / 'no convergence in maxfn iterations in zbrent' /
  data wmess(4)  / 'f(a) and f(b) have the same sign in zbrent' /
  data wmess(5)  / 'negative pivot encountered in factorization' /
  data wmess(6)  / 'unknown warning' /
!
  ier = ierr
  if (ier == 0) return
  if (ier < 0  .and.  level < 0) return
  if (ier > 0  .and.  level < 1) return
  if (ier < -19) ier = -20
  if (ier >   5) ier =   6
  if (ier < 0) write (nout,10)
 10   format (//1x,60('*') /1x,18('*'),' Fatal error ',18('*') /1x,60('*') /)
  if (ier > 0) write (nout,20)
 20   format (//1x,60('*') /1x,16('*'),' Warning error ',16('*') /1x,60('*') /)
  write (nout,23) iname
 23   format (' Routine ',a)
  inum = iabs(ier)
  if (ier > 0) go to 30
!
!  print out fatal errors.
!
  write (nout,25) fmess(inum)
 25   format (1x,a80)
  go to 999
!
!  print out warning errors.
!
 30   write (nout,25) wmess(inum)
  if (inum /= 2) go to 999
  temp = 500.0*srelpr
  write (nout,35) zeta, srelpr, temp
 35   format (1x,'rparm(1) =',e10.3,' (zeta)' &
    / 1x, 'a value this small may hinder convergence' &
    / 1x, 'since machine precision srelpr = ',e10.3 &
    / 1x, 'zeta reset to ',e10.3)
!
!  print ending line.
!
 999  write (nout,1000)
 1000 format (/1x,60('*')/)
  return
end
subroutine fillb (nn,coef,jcoef,iblock,wksp,iwksp,ier)
!
!*******************************************************************************
!
!! FILLB calculates block fill-in for block factorization methods.
!
!     (symmetric diagonal storage)
!
!  Parameters:
!
!       n       order of system
!       coef    real matrix coefficient array
!       jcoef   integer matrix coefficient array
!       iblock  array for block information
!       wksp    real workspace array
!       iwksp   integer workspace array
!       ier     error flag
!
!  
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension coef(1), wksp(1)
  integer   jcoef(1), iblock(3,3), iwksp(1)
!
  n = nn
!
!  determine block fill-in pattern.
!
  if (lvfill > 0) then
     lbhbsa = lbhb
     do 25 lv = 1,lvfill
        lbhbl = lbhb
        do 20 j1 = 3,lbhb
           do 15 j2 = 3,lbhb
              jd = iblock(1,j1) - iblock(1,j2)
              if (jd <= 0) go to 15
              do 10 j3 = 3,lbhbl
                 if (iblock(1,j3) == jd) go to 15
 10               continue
              lbhbl = lbhbl + 1
              iblock(1,lbhbl) = jd
              iblock(3,lbhbl) = 0
 15            continue
 20         continue
        lbhb = lbhbl
 25      continue
  end if
!
!  compute constants and check for sufficient workspace.
!
  call needw ('fillb',1,iblk,3*lbhb,ier)
  if (ier < 0) return
  nwdiag = iblock(3,1)
  nwnew = nwdiag + ltrunc
  iipnt = iblk + 3*lbhb
  ifactr = irpnt
  nwk = 3*lbhb + maxnz + ltrunc + (lbhb-2)*(2*nwnew-1)
  call needw ('fillb',1,iblk,nwk,ier)
  if (ier < 0) return
  do 30 j = 1,nwnew
 30   iwksp(iipnt+j-1) = j - 1
  iblock(3,1) = nwnew
!
!  determine diagonal numbers in filled-in block matrix.
!
  if (lvfill > 0) then
     jmax = 3
     do 32 j = 3,lbhbsa
        if (iblock(1,j) > iblock(1,jmax)) jmax = j
 32      continue
     jnext = iipnt + nwnew
     do 50 jjc = 3,lbhb
        if (jjc <= lbhbsa) then
           jstc = iblock(2,jjc)
           mc = iblock(3,jjc)
           j1 = jnext
           do 35 j = 1,mc
              iwksp(jnext) = jcoef(nwdiag+jstc+j-1)
              jnext = jnext + 1
 35            continue
           j2 = jnext - 1
        end if
        if (jjc == jmax) go to 50
        jblkc = iblock(1,jjc)
        inc = jblkc*kblsz
        lim1 = inc - (nwnew - 1)
        lim2 = inc + (nwnew - 1)
        do 45 j = lim1,lim2
           if (jjc <= lbhbsa) then
              do 40 jj = j1,j2
                 if (iwksp(jj) == j) go to 45
 40               continue
           end if
           iwksp(jnext) = j
           jnext = jnext + 1
           iblock(3,jjc) = iblock(3,jjc) + 1
 45         continue
 50      continue
     if (lbhb >= 4) then
        do 52 jjc = 4,lbhb
 52         iblock(2,jjc) = iblock(2,jjc-1) + iblock(3,jjc-1)
     end if
  end if
!
!  copy matrix into wksp.
!
  if (propa) then
     nfactr = n*nwnew
     nfacti = 3*lbhb
  end if
  if (.not. propa .and. lvfill == 0) then
     nfactr = n*(maxnz + ltrunc)
     nfacti = 3*lbhb
  end if
  if (lvfill > 0) then
     ndg = 0
     do 55 j = 1,lbhb
 55      ndg = ndg + iblock(3,j)
     nfactr = n*ndg
     nfacti = ndg + 3*lbhb
  end if
  call needw ('fillb',0,ifactr,nfactr,ier)
  if (ier < 0) return
  call needw ('fillb',1,ifacti,nfacti,ier)
  if (ier < 0) return
  call vfill (nfactr,wksp(ifactr),0.0)
  ipt1 = 1
  ipt2 = ifactr
  do 60 j = 1,nwdiag
     call vcopy (n,coef(ipt1),wksp(ipt2))
     ipt1 = ipt1 + ndim
     ipt2 = ipt2 + n
 60   continue
  iwkpt2 = ifactr + n*nwnew
  ipt2 = iwkpt2
  if (.not. propa .and. lvfill == 0) then
     do 62 j = nwdiag+1,maxnz
        call vcopy (n,coef(ipt1),wksp(ipt2))
        ipt1 = ipt1 + ndim
        ipt2 = ipt2 + n
 62      continue
  end if
  if (lvfill > 0) then
     j1 = iipnt + nwnew
     j2 = iipnt + ndg - 1
     do 70 j = nwdiag+1,maxnz
        jcol = jcoef(j)
        ipt1 = (j - 1)*ndim + 1
        do 65 jj = j1,j2
           if (iwksp(jj) /= jcol) go to 65
           ipt2 = iwkpt2 + (jj-j1)*n
           call vcopy (n,coef(ipt1),wksp(ipt2))
           go to 70
 65         continue
 70      continue
  end if
  irpnt = ifactr + nfactr
  iipnt = ifacti + nfacti
  return
end
subroutine fillbc (nn,ncolor,coef,jcoef,iblock,wksp,iwksp,ier)
!
!*******************************************************************************
!
!! FILLBC sets up WKSP for block factorization methods.
!
!
!     (multicolor nonsymmetric diagonal)
!
!  Parameters:
!
!       n       order of system
!       coef    real matrix coefficient array
!       jcoef   integer matrix coefficient array
!       iblock  array for block information
!       wksp    real workspace array
!       iwksp   integer workspace array
!       ier     error flag
!
!  
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncol, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension coef(1), wksp(1)
  integer   jcoef(1), iblock(3,ncolor,3), iwksp(1)
!
  n = nn
!
!  compute constants and check for sufficient workspace.
!
  ndt = 0
  ndb = 0
  do 10 j = 1,ncolor
     ndt = max (ndt,iblock(3,j,1)-1)
     ndb = max (ndb,iblock(3,j,2))
 10   continue
  nwdiag = ndt + ndb + 1
  nwnew = nwdiag + 2*ltrunc
  ifactr = irpnt
!
!  copy matrix into wksp.
!
  if (propa) nfactr = n*nwnew
  if (.not. propa) nfactr = n*nwnew + n*(maxd-nwdiag)
  call needw ('fillbc',0,ifactr,nfactr,ier)
  if (ier < 0) return
  call needw ('fillbc',1,iipnt,nwnew*ncolor,ier)
  if (ier < 0) return
  call vfill (nfactr,wksp(ifactr),0.0)
  ipt1 = 1
  ipt2 = ifactr
  do 15 j = 1,ndt+1
     call vcopy (n,coef(ipt1),wksp(ipt2))
     ipt1 = ipt1 + ndim
     ipt2 = ipt2 + n
 15   continue
  ipt2 = ipt2 + n*ltrunc
  do 20 j = ndt+2,nwdiag
     call vcopy (n,coef(ipt1),wksp(ipt2))
     ipt1 = ipt1 + ndim
     ipt2 = ipt2 + n
 20   continue
  iwkpt2 = ifactr + n*nwnew
  ipt2 = iwkpt2
  if (.not. propa) then
     do 25 j = nwdiag+1,maxd
        call vcopy (n,coef(ipt1),wksp(ipt2))
        ipt1 = ipt1 + ndim
        ipt2 = ipt2 + n
 25      continue
  end if
  irpnt = ifactr + nfactr
  do 40 ico = 1,ncolor
     do 30 j = 1,ndt+ltrunc+1
 30      iwksp(iipnt+(j-1)*ncolor+ico-1) = j - 1
     do 35 j = ndt+ltrunc+2,nwnew
 35      iwksp(iipnt+(j-1)*ncolor+ico-1) = -(j - ndt - ltrunc - 1)
 40   continue
  do 45 ico = 1,ncolor
     iblock(3,ico,1) = ndt + ltrunc + 1
     iblock(3,ico,2) = ndb + ltrunc
     iblock(2,ico,2) = iblock(2,ico,1) + iblock(3,ico,1)
 45   continue
  return
end
subroutine fillbn (nn,coef,jcoef,iblock,wksp,iwksp,ier)
!
!*******************************************************************************
!
!! FILLBN calculates block fill-in for block factorization methods.
!
!
!     (nonsymmetric diagonal storage)
!
!  Parameters:
!
!       n       order of system
!       coef    real matrix coefficient array
!       jcoef   integer matrix coefficient array
!       iblock  array for block information
!       wksp    real workspace array
!       iwksp   integer workspace array
!       ier     error flag
!
!  
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension coef(1), wksp(1)
  integer   jcoef(1), iblock(3,3), iwksp(1)
!
  n = nn
!
!  determine block fill-in pattern.
!
  if (lvfill > 0) then
     lbhbsa = lbhb
     do 25 lv = 1,lvfill
        lbhbl = lbhb
        do 20 j1 = 3,lbhb
           do 15 j2 = 3,lbhb
              jd = iblock(1,j1) + iblock(1,j2)
              if (iblock(1,j1)*iblock(1,j2) >= 0) go to 15
              do 10 j3 = 1,lbhbl
                 if (iblock(1,j3) == jd) go to 15
 10               continue
              lbhbl = lbhbl + 1
              iblock(1,lbhbl) = jd
              iblock(3,lbhbl) = 0
 15            continue
 20         continue
        lbhb = lbhbl
 25      continue
  end if
!
!  compute constants and check for sufficient workspace.
!
  call needw ('fillbn',1,iblk,3*lbhb,ier)
  if (ier < 0) return
  ndt = iblock(3,1) - 1
  ndb = iblock(3,2)
  nwdiag = ndt + ndb + 1
  nwnew = nwdiag + 2*ltrunc
  iipnt = iblk + 3*lbhb
  ifactr = irpnt
  nwk = 3*lbhb + maxnz + 2*ltrunc + (lbhb-2)*nwnew
  call needw ('fillbn',1,iblk,nwk,ier)
  if (ier < 0) return
  do 30 j = 1,ndt+ltrunc+1
 30   iwksp(iipnt+j-1) = j - 1
  do 31 j = ndt+ltrunc+2,nwnew
 31   iwksp(iipnt+j-1) = -(j - ndt - ltrunc - 1)
  iblock(3,1) = ndt + ltrunc + 1
  iblock(3,2) = ndb + ltrunc
  iblock(2,2) = iblock(2,1) + iblock(3,1)
!
!  determine diagonal numbers in filled-in block matrix.
!
  if (lvfill > 0) then
     jmax = 3
     jmin = 3
     do 32 j = 3,lbhbsa
        if (iblock(1,j) > iblock(1,jmax)) jmax = j
        if (iblock(1,j) < iblock(1,jmin)) jmin = j
 32      continue
     jnext = iipnt + nwnew
     do 50 jjc = 3,lbhb
        if (jjc <= lbhbsa) then
           jstc = iblock(2,jjc)
           mc = iblock(3,jjc)
           j1 = jnext
           do 35 j = 1,mc
              iwksp(jnext) = jcoef(nwdiag+jstc+j-1)
              jnext = jnext + 1
 35            continue
           j2 = jnext - 1
        end if
        if (jjc == jmax .or. jjc == jmin) go to 50
        jblkc = iblock(1,jjc)
        inc = jblkc*kblsz
        lim1 = inc - (ndb + ltrunc)
        lim2 = inc + (ndt + ltrunc)
        do 45 j = lim1,lim2
           if (jjc <= lbhbsa) then
              do 40 jj = j1,j2
                 if (iwksp(jj) == j) go to 45
 40               continue
           end if
           iwksp(jnext) = j
           jnext = jnext + 1
           iblock(3,jjc) = iblock(3,jjc) + 1
 45         continue
 50      continue
     if (lbhb >= 4) then
        do 52 jjc = 4,lbhb
 52         iblock(2,jjc) = iblock(2,jjc-1) + iblock(3,jjc-1)
     end if
  end if
!
!  copy matrix into wksp.
!
  if (propa) then
     nfactr = n*nwnew
     nfacti = 3*lbhb
  end if
  if (.not. propa .and. lvfill == 0) then
     nfactr = n*(maxnz + 2*ltrunc)
     nfacti = 3*lbhb
  end if
  if (lvfill > 0) then
     ndg = 0
     do 55 j = 1,lbhb
 55      ndg = ndg + iblock(3,j)
     nfactr = n*ndg
     nfacti = ndg + 3*lbhb
  end if
  call needw ('fillbn',0,ifactr,nfactr,ier)
  if (ier < 0) return
  call needw ('fillbn',1,ifacti,nfacti,ier)
  if (ier < 0) return
  call vfill (nfactr,wksp(ifactr),0.0)
  ipt1 = 1
  ipt2 = ifactr
  do 60 j = 1,ndt+1
     call vcopy (n,coef(ipt1),wksp(ipt2))
     ipt1 = ipt1 + ndim
     ipt2 = ipt2 + n
 60   continue
  ipt2 = ipt2 + n*ltrunc
  do 61 j = ndt+2,nwdiag
     call vcopy (n,coef(ipt1),wksp(ipt2))
     ipt1 = ipt1 + ndim
     ipt2 = ipt2 + n
 61   continue
  iwkpt2 = ifactr + n*nwnew
  ipt2 = iwkpt2
  if (.not. propa .and. lvfill == 0) then
     do 62 j = nwdiag+1,maxnz
        call vcopy (n,coef(ipt1),wksp(ipt2))
        ipt1 = ipt1 + ndim
        ipt2 = ipt2 + n
 62      continue
  end if
  if (lvfill > 0) then
     j1 = iipnt + nwnew
     j2 = iipnt + ndg - 1
     do 70 j = nwdiag+1,maxnz
        jcol = jcoef(j)
        ipt1 = (j - 1)*ndim + 1
        do 65 jj = j1,j2
           if (iwksp(jj) /= jcol) go to 65
           ipt2 = iwkpt2 + (jj-j1)*n
           call vcopy (n,coef(ipt1),wksp(ipt2))
           go to 70
 65         continue
 70      continue
  end if
  irpnt = ifactr + nfactr
  iipnt = ifacti + nfacti
  return
end
subroutine filln (maxnz,jcoef)
!
!*******************************************************************************
!
!! FILLN determines the fill-in diagonals for nonsymmetric diagonal storage.
!
!
!  Parameters:
!
!        maxnz   upon input, the number of diagonals
!                upon output, the number of diagonals with fill-in
!        jcoef   upon input, the diagonal numbers
!                upon output, the diagonal numbers with fill-in
!
!  
!
  integer jcoef(2)
!
  maxn = maxnz
  do 20 j1 = 1,maxnz
     do 15 j2 = 1,maxnz
        jd = jcoef(j1) + jcoef(j2)
        if (jcoef(j1)*jcoef(j2) >= 0) go to 15
        do 10 j3 = 1,maxn
           if (jcoef(j3) == jd) go to 15
 10         continue
        maxn = maxn + 1
        jcoef(maxn) = jd
 15      continue
 20   continue
  maxnz = maxn
  return
end
subroutine fillnp (ndim,nn,maxcc,jc,c,mwidth,ier)
!
!*******************************************************************************
!
!! FILLNP determines the fill-in structure.
!
!
!     (Purdue storage, nonsymmetric matrix)
!
!  Parameters:
!
!          ndim   row dimension of jc and c arrays
!          n      order of system (= nn)
!          maxc   upon input, maxc is the number of columns in
!                  the c array
!                 upon output, maxc is the number of columns in
!                  the c array with fill-in
!          jc     integer array of active size n by maxc giving the
!                  column numbers of the corresponding elements in c
!          c      array of active size n by maxc giving the
!                  coefficients of the off-diagonal elements
!          mwidth maximum column width to be allowed for fill-in
!          ier    error code
!                  =    0  no errors detected
!                  =   -2  mwidth too small to accomodate fill-in
!
!  
!
  integer   jc(ndim,1)
  dimension c(ndim,1)
!
!
  n = nn
  maxc = maxcc
  maxu = maxc
!
  if (maxc < 1) return
  nm1 = n - 1
  do 45 k = 1,nm1
     kp1 = k + 1
     do 40 j1 = 1,maxc
     do 35 i = kp1,n
        if (jc(i,j1) /= k) go to 35
        do 30 j2 = 1,maxc
           j = jc(k,j2)
           if (j <= k .or. j == i) go to 30
           do 10 j3 = 1,maxu
              if (j == iabs(jc(i,j3))) go to 30
 10            continue
           do 15 j3 = 1,maxu
              if (jc(i,j3) /= i) go to 15
              jc(i,j3) = -j
              go to 30
 15            continue
           maxu = maxu + 1
           if (maxu <= mwidth) go to 20
           ier = -2
           return
 20            do 25 ii = 1,n
              jc(ii,maxu) = ii
              c(ii,maxu)= 0.0
 25            continue
           jc(i,maxu) = -j
 30         continue
 35      continue
 40      continue
 45   continue
!
!  decode new elements of jt, jb.
!
  do 55 j = 1,maxu
     do 50 i = 1,n
 50      jc(i,j) = iabs(jc(i,j))
 55   continue
  maxcc = maxu
  return
end
subroutine fills (maxt,jt)
!
!*******************************************************************************
!
!! FILLS determines the fill-in diagonals for symmetric diagonal storage.
!
!
!  Parameters:
!
!        maxt    upon input, the number of diagonals in the
!                 upper triangle
!                upon output, the number of diagonals in the
!                 upper triangle with fill-in
!        jt      upon input, the diagonal numbers in the upper
!                 triangle
!                upon output, the diagonal numbers in the upper
!                 triangle with fill-in
!
!  
!
  integer jt(1)
!
  maxn = maxt
  do 20 j1 = 1,maxt
     do 15 j2 = 1,maxt
        jd = jt(j1) - jt(j2)
        if (jd <= 0) go to 15
        do 10 j3 = 1,maxn
           if (jt(j3) == jd) go to 15
 10         continue
        maxn = maxn + 1
        jt(maxn) = jd
 15      continue
 20   continue
  maxt = maxn
  return
end
subroutine fillsp (ndim,nn,maxtt,jt,t,mwidth,ier)
!
!*******************************************************************************
!
!! FILLSP determines the fill-in structure.
!
!     (Purdue storage, symmetric matrix)
!
!  Parameters:
!
!          ndim   row dimension of t and jt arrays
!          n      order of system (= nn)
!          maxt   upon input, maxt is the number of columns in
!                  the t array
!                 upon output, maxt is the number of columns in
!                  the t array with fill-in
!          jt     integer array of active size n by maxt giving the
!                  column numbers of the corresponding elements in t
!          t      array of active size n by maxt giving the
!                  coefficients of the upper triangle of the matrix
!          mwidth maximum column width of jt and t to be allowed
!          ier    error code
!                  =   0     no error detected
!                  =  -2     mwidth too small to store factor
!
!  
!
  dimension t(ndim,1)
  integer   jt(ndim,1)
!
!
  n = nn
  maxt = maxtt
  maxu = maxt
  ier = 0
!
  if (maxt < 1) return
  nm1 = n - 1
  do 40 k = 1,nm1
     do 35 j1 = 1,maxt
        jcol1 = jt(k,j1)
        if (jcol1 <= 0 .or. jcol1 == k) go to 35
        do 30 j2 = 1,maxt
           jcol2 = jt(k,j2)
           if (jcol2 <= 0 .or. jcol2 == k) go to 30
           if (jcol2 <= jcol1) go to 30
           do 10 j3 = 1,maxu
              if (jcol2 == iabs(jt(jcol1,j3))) go to 30
 10            continue
           do 15 j3 = 1,maxu
              if (jt(jcol1,j3) /= jcol1) go to 15
              jt(jcol1,j3) = -jcol2
              go to 30
 15            continue
           maxu = maxu + 1
           if (maxu <= mwidth) go to 20
           ier = -2
           return
 20            do 25 i = 1,n
              jt(i,maxu) = i
              t(i,maxu) = 0.0
 25            continue
           jt(jcol1,maxu) = -jcol2
 30         continue
 35      continue
 40   continue
!
!  decode new elements of jt.
!
  do 50 j = 1,maxu
     do 45 i = 1,n
 45      jt(i,j) = iabs(jt(i,j))
 50   continue
  maxtt = maxu
  return
end
subroutine gauss (ndim,n,a,rhs,u,ier)
!
!*******************************************************************************
!
!! GAUSS is a Gaussian elimination routine.
!
  dimension a(ndim,ndim), rhs(ndim), u(ndim)
  common / itcom4 / keygs, srelpr, keyzer
  ier = 0
  if (n == 1) go to 190
  do 1 i = 1,n-1
  if (abs(a(i,i)) < srelpr**2) go to 999
  do 10 j = i+1,n
  fact = a(j,i)/a(i,i)
  a(j,i) = 0.e0
  do 2 k = i+1,n
 2    a(j,k) = a(j,k) - fact*a(i,k)
  rhs(j) = rhs(j) - fact*rhs(i)
 10   continue
 1    continue
!
 190  do 3 i = 1,n
  k = n - i + 1
  if (abs(a(k,k)) < srelpr**2) go to 999
  u(k) = rhs(k)
  if (i == 1) go to 44
  do 4 j = k+1,n
  u(k) = u(k) - u(j)*a(k,j)
 4    continue
 44   u(k) = u(k)/a(k,k)
 3    continue
  return
 999  ier = -100
  return
end
subroutine getblk (coef,jcoef,n,nblk,nband,ctac,nw,ier)
!
!*******************************************************************************
!
!! GETBLK computes and factors the matrix (C**t)*A*C and factors it.
!
!
!  this is a utility routine for the cgcr algorithm.
!  here, each column of c is zero
! everywhere except it is all 1's on one of its blocks.
!
  dimension ctac(nblk,1), coef(1), jcoef(2)
  logical symm
!
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
!
!
  nband = 0
!
!  find the bandwidth
!
!
  idmin = 0
  idmax = 0
  do 10 i=1,maxnz
  idiag = jcoef(i)
  idmin = min (idmin,idiag)
 10   idmax = max (idmax,idiag)
  if (nstore == 2) idmin = - idmax
  ihalf = max(-idmin,idmax)
  nbsiz = n / nblk
  nhband = (ihalf+nbsiz-1)/nbsiz
  nband = 1 + 2*nhband
!
! now form the matrix.  basically what we need to do here is to
! add up all the elements in each block of the a-matrix.
!
  if (nblk*nband > nw) go to 999
  nw = nblk*nband
!
  call vfill (nblk*nband,ctac,0e0)
!
! loop over the diagonals.
!
  do 1 i=1,maxnz
  idiag = jcoef(i)
  ibeg = max (1,1-idiag)
  iend = min (n,n-idiag)
  ibbeg = 1 + (ibeg-1)/nbsiz
  ibend = 1 + (iend-1)/nbsiz
  ibase = (i-1)*ndim
!
  symm = nstore==2 .and. idiag/=0
  idm1 = idiag - 1
  iomid = -idm1
  nmid = n - idiag
  nhbp1 = nhband + 1
! loop over the rows of ctac.
!
  do 2 j=ibbeg,ibend
  ibeg = max(1+(j-1)*nbsiz,iomid)
  iend = min(j*nbsiz,nmid)
!     ic1 = (ibeg+idiag-1)/nbsiz + 1
!     ic2 = (iend+idiag-1)/nbsiz + 1
!     id1 = ic1 - j + nhband + 1
!     id2 = ic2 - j + nhband + 1
  itemp1 = (ibeg+idm1)/nbsiz
  itemp2 = (iend+idm1)/nbsiz
  id1 = itemp1 + 2 - j + nhband
  id2 = itemp2 + 2 - j + nhband
  j1s = j + id1 - nhbp1
  j2s = j + id2 - nhbp1
  id1s = nband - id1 + 1
  id2s = nband - id2 + 1
  if (id1 /= id2) go to 3
!     ctac(j,id1) = ctac(j,id1)
!    a              + vadd(iend-ibeg+1,coef(ibase+ibeg))
  do 41 ii=ibeg,iend
  if (symm) ctac(j1s,id1s) = ctac(j1s,id1s) + coef(ibase+ii)
 41   ctac(j,id1) = ctac(j,id1) + coef(ibase+ii)
  go to 2
!3    imid = 1 + (ic2-1)*nbsiz - idiag
 3    imid = iomid + itemp2*nbsiz
!     ctac(j,id1) = ctac(j,id1)
!    a              + vadd(imid-ibeg  ,coef(ibase+ibeg))
  do 42 ii=ibeg,imid-1
  if (symm) ctac(j1s,id1s) = ctac(j1s,id1s) + coef(ibase+ii)
 42   ctac(j,id1) = ctac(j,id1) + coef(ibase+ii)
!     ctac(j,id2) = ctac(j,id2)
!    a              + vadd(iend-imid+1,coef(ibase+imid))
  do 43 ii=imid,iend
  if (symm) ctac(j2s,id2s) = ctac(j2s,id2s) + coef(ibase+ii)
 43   ctac(j,id2) = ctac(j,id2) + coef(ibase+ii)
!
 2    continue
!
 1    continue
!
!  do lu factorization
!
  do 31 i=1,nblk-1
  denom = ctac(i,nhbp1)
  if (abs(denom) < srelpr) go to 998
  xpivot = 1e0 / denom
  nsubmt = min(nhband,nblk-i)
  do 30 j=1,nsubmt
  ipj = i + j
  ind2 = nhbp1 - j
  do 30 k=1,nsubmt
!30   ctac(i+j,nhband-j+1+k) = ctac(i+j,nhband-j+1+k)
!    a  - xpivot*ctac(i+j,nhband-j+1)*ctac(i,nhband+1+k)
  ind = nhbp1 - j + k
 30   ctac(ipj,ind) = ctac(ipj,ind) - xpivot*ctac(ipj,ind2)*ctac(i,nhbp1+k)
  do 32 j=1,nsubmt
  ipj = i + j
  ind1 = nhbp1 - j
  ind2 = nhbp1 + j
!     ctac(i+j,nhband+1-j) = ctac(i+j,nhband+1-j)*xpivot
!32   ctac(i  ,nhband+1+j) = ctac(i  ,nhband+1+j)*xpivot
  ctac(ipj,ind1) = ctac(ipj,ind1)*xpivot
 32   ctac(i  ,ind2) = ctac(i  ,ind2)*xpivot
 31   continue
  return
!
!
! error returns
!
! breakdown.
!
 998  ier = -6
  call ershow (ier,'getblk')
  return
!
! insuff. memory.
!
 999  ier = -2
  call ershow (ier,'getblk')
  nw = nblk*nband
  return
end
subroutine gmres (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef, &
  n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! GMRES is the user interface to the truncated/restarted GMRES algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2)
  dimension wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call gmresw (suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs, &
    wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine gmresw (suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
  wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! GMRESW runs the truncated/restarted GMRES algorithm.  
!
!
!  a detailed
! description of this useful algorithm may be found in the paper,
! "gmres: a generalized minimal residual algorithm for solving
! nonsymmetric linear systems", youcef saad and martin h. schultz,
! siam j. sci. stat. comput., v. 7, no. 3, july 1986.
!
! further scoop on how to set up qr factorizations can be obtained in
! "practical use of some krylov subspace methods for solving
! indefinite and unsymmetric linear systems", youcef saad, siam j. sci.
! stat. comput., v. 5, no. 1, march 1984.
!
! the advantage of this algorithm over its competitors orthomin and gcr
! is that work and storage are saved by avoiding the computation of
! certain vectors.
!
! this routine now handles right and 2-sided preconditioning.  the main
! thing to note about this is that a new table of basis vecttors is now
! necessary, to use to update the solution.
!
! this routine also avoids explicit scaling of the p and w vectors.
!
! for the pure restarted case, we actually compute the final arnoldi
! vector, rather than just estimating its norm.  this is a diversion
! from the saad/schultz paper.  this was done because in some cases it
! was found that the norm estimation was subject to significant
! numerical error.
!
! modified feb. 1990 to make the restarted method more efficient.
! specifically, new formulas were installed for the scalar part of
! the computation to give an optimal asymptotic dependence on ns2.
!
  dimension coef(*), jcoef(*), wfac(*), jwfac(*)
  dimension u(*), ubar(*), rhs(*), wk(*)
  logical uneed, zneed
  external suba, subql, subqr
  dimension iparm(30), rparm(30)
  logical iql, iqr
  logical trunc, exact, rstrt, rstrtd, zhvold
  logical havest, hadest, evadpt
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
!  indexing functions.
!
! the following function accesses the arnoldi vectors.
  indp(i) = ip + mod(i,nv)*n
!
! the following accesses q-r times the arnoldi vectors
  indpt(i) = ipt + mod(i,nvt)*n
!
! fudge factor for the arnoldi vectors.  p(actual) = p(stored)*pfudge.
! (we do the same trick with a*p.)
  indpf(i) = ipf + mod(i,nv)
!
! the following accesses the w-vectors.
  indw(i) = iw + n*mod(i,nv)
!
! fudge factors for the w vectors.
!
! (similarly, the vector "xi" is fudged.)
  indwf(i) = iwf + mod(i,nv)
!
! the following accesses the Hessenberg matrix -- stored by diagonals.
!
  indhes(i,j) = ihess + (i-1) + (j-i+1)*nhess
!
! the following are the cosines and sines of the rotations.
  indc(i) = icos + mod(i,nrot)
  inds(i) = isin + mod(i,nrot)
!
! the following accesses the u matrix -- stored by columns.
!
  indu(i,j) = iu + j-i+1 + mod(j-1,nuc)*nbwuh
!
! the following accesses the z-vector.
!
  indzc(i) = izc + mod(i-1,nzc)
!
!  preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 11
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 996
  iql = iqlr==1 .or. iqlr==3
  iqr = iqlr==2 .or. iqlr==3
  iadpt = ns3
  evadpt = (maxadd.or.minadd) .and. iadpt/=0
  trunc = ns1 < (ns2-1)
  exact = .not. trunc
  if (ns1 < 2) go to 995
  if (level >= 2) write (nout,496)
496   format (' gmres')
!
!  initialize the stopping test.
!
  call inithv (0)
  zdhav = .not. (trunc .and. .not.exact)
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
! uneed tells us whether u must be computed explicitly per iteration.
! similarly for zneed.
  uneed = rcalp .or. udhav .or. ntest == 6 .or. level >= 3
  zneed = zcalp
  hadest = .false.
!
!  associated integer variables.
!
!  effective ns2.
  ns2e = min(ns2,itmax)
!  length of diags of hess matrix.
  nhess = ns2e + 2
!  bandwidth of the hess matrix.
  nbwh = min(ns1+1,ns2e+1)
!  bandwidth of u-or-h.
  nbwuh = min(ns1+2, ns2e+1)
!  number columns stored of the u matrix.
  if (     trunc) nuc = 1
  if (.not.trunc) nuc = ns2e
!  size of arnoldi-vector tables.
  nv = min(ns1,ns2e+1)
  nvt = nv
  if (iqr .and. .not.trunc) nvt = nv - 1
  if (iqr .and.      trunc) nvt = 1
!  number of givens rotations to store.
  nrot = min(ns1,ns2e)
!  number of elts of z-vector to store.
  if (     trunc) nzc = 2
  if (.not.trunc) nzc = ns2e + 1
!
!  memory layout.
!
  ihess = 1
                   ipt = ihess + nhess*nbwh
  if (.not.evadpt) ipt = ihess
                ip = ipt + n*nvt
  if (.not.iqr) ip = ipt
  izc = ip + n*nv
  icos = izc + nzc
  isin = icos + nrot
  iy = isin + nrot
                             iu = iy + ns2e
  if (trunc .or. .not.uneed) iu = iy
  ipz = iu + nbwuh*nuc
             ipf = ipz + ns2e+1
  if (trunc) ipf = ipz
  iz = ipf + nv
  iw = iz + n
                   iwf = iw + n*nv
                   ixi = iwf + nv
                   iv1 = ixi + n
  if (.not. trunc) iv1 = iw
  iv2 = iv1 + n
  nwusd = max(nwusd,iv2+n-1)
!
! check the memory usage.
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  rstrtd = .true.
!
!  begin iteration loop.
!
! handle first iteration after restart.
!
!
 10   call inithv (1)
  zdhav = .not.(trunc.and..not.exact) .and. in/=0
  if (.not. rstrtd) go to 100
!  get resid.
  if (.not. zhave) then
    if (iql) then
      call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
      call vexopy (n,wk(iv1),rhs,wk(iv1),2)
      call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iz))
    else
      call suba (coef,jcoef,wfac,jwfac,n,u,wk(iz))
      call vexopy (n,wk(iz),rhs,wk(iz),2)
    end if
    zhave = .true.
  end if
!  get resid norm.
  if (.not. zdhav) then
    zdot = vdot (n,wk(iz),wk(iz))
    zdhav = .true.
  end if
  if (zdot < 0e0) go to 994
  vnorm = sqrt(zdot)
  if (vnorm < srelpr**2) go to 997
  call vcopy (n,wk(iz),wk(indp(is)))
  wk(indpf(is)) = 1e0/vnorm
  wk(indzc(is+1)) = vnorm
!
!  perform stopping test.
!
 100  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,wk(iz),xxx,wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
! rstrt tells us whether this is the last step before restarting.
!
  rstrt = (is+1 == ns2)
  if (evadpt .and. is==0) call vfill (nhess*nbwh,wk(ihess),0e0)
!
!  compute the new arnoldi vector.
!
! pn(is+1)*p(is+1) = a*p(is) + sum (i=0 to is) (beta(is+1,i)*p(i)),
!
!  get a times old vec.
  if (iqr) call subqr (coef,jcoef,wfac,jwfac,n,wk(indp (is)),wk(indpt(is)))
  if (iql) then
    call suba (coef,jcoef,wfac,jwfac,n,wk(indpt(is)),wk(iv1))
    call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  else
    call suba (coef,jcoef,wfac,jwfac,n,wk(indpt(is)),wk(iv2))
  end if
  apf = wk(indpf(is))
!  compute arnoldi vector.
  ibeg = max(is+1-ns1,0)
  iend = is
  if (ibeg > 0) wk(indu(ibeg,is+1)) = 0e0
  pfnew = apf
  do 199 i = ibeg,iend
    h = vdot (n,wk(indp(i)),wk(iv2)) * wk(indpf(i))*apf
    wk(indu(i+1,is+1)) = h
    if (evadpt) wk(indhes(i+1,is+1)) = h
    if (i == ibeg) call vtriad (n,wk(indp(is+1)),wk(iv2), &
                                -h*wk(indpf(i))/pfnew,wk(indp(i)),1)
    if (i /= ibeg) call vtriad (n,wk(indp(is+1)),wk(indp(is+1)), &
                                -h*wk(indpf(i))/pfnew,wk(indp(i)),1)
 199  continue
  wk(indpf(is+1)) = pfnew
!  get norm.
  dot = vdot (n,wk(indp(is+1)),wk(indp(is+1))) * pfnew**2
  vnorm = sqrt(dot)
  if (vnorm < srelpr**2) go to 192
  wk(indu(is+2,is+1)) = vnorm
  if (evadpt) wk(indhes(is+2,is+1)) = vnorm
!  scale.
  wk(indpf(is+1)) = wk(indpf(is+1))/vnorm
  if (abs(wk(indpf(is+1)))<srelpr .or. &
    abs(wk(indpf(is+1)))>1e0/srelpr) then
    call vtriad (n,wk(indp(is+1)),xxx,wk(indpf(is+1)), wk(indp(is+1)),2)
    wk(indpf(is+1)) = 1e0
  end if
!
!  update the qr factorization.
!
 192  continue
!  apply old rotations.
  ibgn = max(0,is-ns1)
  iuold = indu(ibgn+1,is+1)
  do 7977 i = ibgn, is-1
  iunew = indu(i+2,is+1)
  ut = wk(iuold)
  h  = wk(iunew)
  ctmp = wk(indc(i+1))
  stmp = wk(inds(i+1))
  wk(iuold) =  ctmp*ut + stmp*h
  wk(iunew) = -stmp*ut + ctmp*h
 7977 iuold = iunew
  iunew = indu(is+2,is+1)
!  calc new rotation.
  v1 = wk(iuold)
  v2 = wk(iunew)
  denom = sqrt (v1**2 + v2**2)
  if (denom < srelpr) go to 998
  wk(indc(is+1)) = v1/denom
  wk(inds(is+1)) = v2/denom
!  apply new rotation.
  wk(iuold) = denom
  wk(iunew) = 0e0
!
!  compute w, if needed.
!
  uc = wk(indu(is+1,is+1))
  if (abs(uc) < srelpr**2) go to 998
  if (.not.trunc) go to 572
!
!  case of explicit w calc.
!
  if (is == 0) then
    call vcopy (n,wk(indpt(is)),wk(indw(1)))
    wk(indwf(is+1)) = wk(indpf(is))/uc
!       call vtriad (n,wk(indw(is+1)),xxx,1.0/uc,wk(indpt(is)),2)
    go to 572
  end if
  wfnew = wk(indpf(is))
  ibeg = max(1,is+1-ns1)
  iend = is
  do 574 i = ibeg, iend
  if (i == ibeg) call vtriad (n,wk(indw(is+1)),wk(indpt(is)), &
       -wk(indu(i,is+1))*wk(indwf(i))/wfnew,wk(indw(i)),1)
  if (i /= ibeg) call vtriad (n,wk(indw(is+1)),wk(indw(is+1)), &
     -wk(indu(i,is+1))*wk(indwf(i))/wfnew,wk(indw(i)),1)
 574  continue
  wk(indwf(is+1)) = wfnew/uc
  if (abs(wk(indwf(is+1)))<srelpr .or.abs(wk(indwf(is+1)))>1e0/srelpr) then
    call vtriad (n,wk(indw(is+1)),xxx,wk(indwf(is+1)),wk(indw(is+1)),2)
    wk(indwf(is+1)) = 1e0
  end if
 572  continue
!
!  get new zc entries.
  wk(indzc(is+2)) = -wk(inds(is+1))*wk(indzc(is+1))
  wk(indzc(is+1)) =  wk(indc(is+1))*wk(indzc(is+1))
!
!  u-vector computation section.
!
  if (trunc) then
!
!  truncated case.
    call vtriad (n,u,u,wk(indzc(is+1))*wk(indwf(is+1)),wk(indw(is+1)),1)
  else
!
!  non-truncated case.
    if (.not.(uneed .or. rstrt)) go to 410
    iynew = iv1
    nwusd = max(nwusd,iynew+ns2e-1)
    if (nwusd > nw) go to 999
!  do back solve on u-matrix.
    nm = is + 1
    do 623 i = nm, 1, -1
    sum = wk(indzc(i))
    do 624 j = i+1, nm
 624    sum = sum - wk(iynew-1+j)*wk(indu(i,j))
 623    wk(iynew-1+i) = sum/wk(indu(i,i))
!  form iterate.
    do 625 i = 0, nm-1
    val = wk(iynew+i)
    if (uneed .and. i/=nm-1) val = val - wk(iy+i)
 625    call vtriad (n,u,u,val*wk(indpf(i)),wk(indpt(i)),1)
    if (uneed) call vcopy (nm,wk(iynew),wk(iy))
  end if
 410  continue
!
!  residual computation section.
!
  zhvold = zhave
  zhave = .false.
  if (trunc) go to 671
!
!  non-truncated case.
!
! do it if resid needed by pstop or if restarting.
  if (zneed .or. rstrt) then
    ipznew = iv1
    nwusd = max(nwusd,ipznew+ns2e)
    if (nwusd > nw) go to 999
    call vcopy (is+1,wk(izc),wk(ipznew))
    wk(ipznew+is+1) = 0e0
!  apply rotations.
    do 644 i = is+1, 1, -1
      v1 = wk(indc(i))*wk(ipznew+i-1) - wk(inds(i))*wk(ipznew+i)
      v2 = wk(inds(i))*wk(ipznew+i-1) + wk(indc(i))*wk(ipznew+i)
      wk(ipznew+i-1) = v1
 644      wk(ipznew+i)   = v2
!  form resid.
    do 645 i = 0, is+1
    val = wk(ipznew+i)
    if (zhvold .and. i/=is+1) val = val - wk(ipz)
 645    call vtriad (n,wk(iz),wk(iz),-val*wk(indpf(i)),wk(indp(i)),1)
    call vcopy (is+2,wk(ipznew),wk(ipz))
    zhave = .true.
  end if
  go to 425
!
!  truncated case.
!
! do it if pstop needs it or if we may restart later.
 671  if ( zneed .or. (itmax>ns2) ) then
!  update xi.
    if (is == 0) then
      call vcopy (n,wk(indp(is)),wk(ixi))
      xif = wk(indpf(is))
    else
      xif = xif*(-wk(inds(is)))
      call vtriad (n,wk(ixi),wk(ixi),wk(indc(is))*wk(indpf(is))/xif, &
        wk(indp(is)),1)
    end if
    if (abs(xif)<srelpr .or. abs(xif)>1e0/srelpr) then
      call vtriad (n,wk(ixi),xxx,xif,wk(ixi),2)
      xif = 1e0
    end if
!  form resid.
    call vtriad (n,wk(iz),wk(iz), &
      -wk(indzc(is+1))*wk(indc(is+1))*xif,wk(ixi),1)
    call vtriad (n,wk(iz),wk(iz), &
      -wk(indzc(is+1))*wk(inds(is+1))*wk(indpf(is+1)),wk(indp(is+1)),1)
    zhave = .true.
  end if
 425  continue
!
!  get resid norm.
!
  if (exact) then
    zdot = wk(indzc(is+2))**2
  end if
!
!  ev est.
!
  if (evadpt) then
    nwhe = nw - (iv1-1)
    call hesest (wk(ihess),nhess,nv+2,is+1,iadpt,havest,emaxnw,eminnw, &
      wk(iv1),nwhe,ier)
    nwusd = max(nwusd,iv1-1+nwhe)
    if (ier /= 0) go to 996
    if (.not. havest) go to 874
    if (hadest) go to 876
    if (maxadd) emax = emaxnw
    if (minadd) emin = eminnw
    hadest = .true.
    go to 874
 876    if (maxadd) emax = amax1 (emax,emaxnw)
    if (minadd) emin = amin1 (emin,eminnw)
  end if
!
!  finish up the iteration.
!
 874  in = in + 1
  is = is + 1
  if (rstrt) is = 0
  rstrtd = rstrt
  go to 10
!
!  wrap it up.
!
!  form u, if not up-to-date.
!
 900  if (uneed .or. rstrtd .or. trunc) go to 901
    iynew = iv1
    nwusd = max(nwusd,iynew+ns2e-1)
    if (nwusd > nw) go to 999
!  do back solve on u-matrix.
    nm = is
    do 663 i = nm, 1, -1
    sum = wk(indzc(i))
    do 664 j = i+1, nm
 664    sum = sum - wk(iynew-1+j)*wk(indu(i,j))
 663    wk(iynew-1+i) = sum/wk(indu(i,i))
!  form iterate.
    do 665 i = 0, nm-1
    val = wk(iynew+i)
 665    call vtriad (n,u,u,val*wk(indpf(i)),wk(indpt(i)),1)
!
!  Head out of here.
!
 901  continue
  if (halt) go to 715
  ier = 1
  call ershow (ier,'gmresw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' gmres converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
!  Error returns.
!
 994  ier = -15
  call ershow (ier,'gmresw')
  return
!
 995  ier = -16
  call ershow (ier,'gmresw')
  return
!
 996  call ershow (ier,'gmresw')
  go to 735
!
 997  ier = -13
  call ershow (ier,'gmresw')
  go to 725
!
 998  ier = -14
  call ershow (ier,'gmresw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'gmresw')
  go to 735
!
end
subroutine hesest (hess,nhess,nd,esize,imode,havest,emax,emin,wk,nw,ier)
!
!*******************************************************************************
!
!! HESEST calculates the extremal eigenvalue moduli of a banded Hessenberg matrix.
!
!
! hess      - the Hessenberg matrix, stored by diagonals
! nhess, nd - dimensions of array hess
! esize     - indicator of how many rows/cols of hess have been
!             filled out so far
! imode     - style of eigenvalue estimation:
!   abs(imode)  - use this size of principal submatrix to do estimate
!   sign(imode) - use either leading or trailing principal submatrix
!
!
  dimension hess(nhess,nd), wk(1)
  logical havest
  integer esize
!
  havest = .false.
  if (imode > 0 .and. esize > imode) return
!
! memory allocation
!
  ndim = min(esize,iabs(imode))
  if (ndim <= 0) return
  imat = 1
  ireal = imat + ndim*ndim
  iimag = ireal + ndim
!
  nwusd = iimag - 1 + ndim
  if (nwusd > nw) go to 999
  nw = nwusd
!
! make the hess matrix into a full matrix
!
  if (imode < 0) go to 1
  ibeg = 1
  iend = esize
  go to 2
 1    ibeg = max (1,esize-iabs(imode)+1)
  iend = esize
 2    call vfill (ndim*ndim,wk(imat),0e0)
  do 3 i=ibeg,iend
  jbeg = max (ibeg,i-1)
  jend = min (ibeg-1+ndim,i+nd-2)
  do 3 j=jbeg,jend
 3    wk(imat+(i-ibeg)+(j-ibeg)*ndim) = hess(i,j-i+2)
!
! call to eispack to calculate eigenvalues
!
  ierr = 0
  call hqr (ndim,ndim,1,ndim,wk(imat),wk(ireal),wk(iimag),ierr)
  if (ierr /= 0) go to 998
!
! find eigenvalues with largest and smallest modulus
!
  emax = wk(ireal)**2 + wk(iimag)**2
  emin = emax
  if (ndim == 1) go to 5
  do 6 i=2,ndim
  vmod = wk(ireal-1+i)**2 + wk(iimag-1+i)**2
  emax = amax1 (emax,vmod)
 6    emin = amin1 (emin,vmod)
!
 5    emax = sqrt (emax)
  emin = sqrt (emin)
  havest = .true.
  return
!
!
! error returns
!
! error in call to eispack
 998  ier = -18
  call ershow (ier,'hesest')
  return
!
! insuff. real workspace
 999  ier = -2
  nw = nwusd
  call ershow (ier,'hesest')
  return
end
subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
!
!*******************************************************************************
!
!! HQR finds the eigenvalues of a real upper Hessenberg matrix by the QR method.
!
!
!  HQR is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          routine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper Hessenberg matrix.  information about
!          the transformations used in the reduction to Hessenberg
!          form by  elmhes  or  orthes, if performed, is stored
!          in the remaining triangle under the Hessenberg matrix.
!
!     on output
!
!        h has been destroyed.  therefore, it must be saved
!          before calling  hqr  if subsequent calculation and
!          back transformation of eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
  integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
  real h(nm,n),wr(n),wi(n)
  real p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
  logical notlas
!
  ierr = 0
  norm = 0.0
  k = 1
!
!  Store roots isolated by balanc and compute matrix norm.
!
  do 50 i = 1, n

     do 40 j = k, n
   40    norm = norm + abs(h(i,j))

     k = i
     if (i >= low .and. i <= igh) go to 50
     wr(i) = h(i,i)
     wi(i) = 0.0
   50 continue

  en = igh
  t = 0.0
  itn = 30*n
!
!  Search for next eigenvalues.
!
   60 if (en < low) go to 1001
  its = 0
  na = en - 1
  enm2 = na - 1
!
!  Look for single small sub-diagonal element
!                for l=en step -1 until low do --........
   70 do 80 ll = low, en
     l = en + low - ll
     if (l == low) go to 100
     s = abs(h(l-1,l-1)) + abs(h(l,l))
     if (s == 0.0) s = norm
     tst1 = s
     tst2 = tst1 + abs(h(l,l-1))
     if (tst2 == tst1) go to 100
   80 continue
!    ........ form shift.
  100 x = h(en,en)
  if (l == en) go to 270
  y = h(na,na)
  w = h(en,na) * h(na,en)
  if (l == na) go to 280
  if (itn == 0) go to 1000
  if (its /= 10 .and. its /= 20) go to 130
!    ........ form exceptional shift.
  t = t + x
!
  do 120 i = low, en
  120 h(i,i) = h(i,i) - x
!
  s = abs(h(en,na)) + abs(h(na,enm2))
  x = 0.75 * s
  y = x
  w = -0.4375 * s * s
  130 its = its + 1
  itn = itn - 1
!    ........ look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do --........
  do 140 mm = l, enm2
     m = enm2 + l - mm
     zz = h(m,m)
     r = x - zz
     s = y - zz
     p = (r * s - w) / h(m+1,m) + h(m,m+1)
     q = h(m+1,m+1) - zz - r - s
     r = h(m+2,m+1)
     s = abs(p) + abs(q) + abs(r)
     p = p / s
     q = q / s
     r = r / s
     if (m == l) go to 150
     tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
     tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
     if (tst2 == tst1) go to 150
  140 continue
!
  150 mp2 = m + 2
!
  do 160 i = mp2, en
     h(i,i-2) = 0.0
     if (i == mp2) go to 160
     h(i,i-3) = 0.0
  160 continue
!    ........ double qr step involving rows l to en and
!                columns m to en........
  do 260 k = m, na
     notlas = k /= na
     if (k == m) go to 170
     p = h(k,k-1)
     q = h(k+1,k-1)
     r = 0.0
     if (notlas) r = h(k+2,k-1)
     x = abs(p) + abs(q) + abs(r)
     if (x == 0.0) go to 260
     p = p / x
     q = q / x
     r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
     if (k == m) go to 180
     h(k,k-1) = -s * x
     go to 190
  180    if (l /= m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
     x = p / s
     y = q / s
     zz = r / s
     q = q / p
     r = r / p
     if (notlas) go to 225
!    ........ row modification.
     do 200 j = k, n
        p = h(k,j) + q * h(k+1,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
  200    continue
!
     j = min(en,k+3)
!    ........ column modification.
     do 210 i = 1, j
        p = x * h(i,k) + y * h(i,k+1)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
  210    continue
     go to 255
  225    continue
!    ........ row modification.
     do 230 j = k, n
        p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
        h(k,j) = h(k,j) - p * x
        h(k+1,j) = h(k+1,j) - p * y
        h(k+2,j) = h(k+2,j) - p * zz
  230    continue
!
     j = min(en,k+3)
!    ........ column modification.
     do 240 i = 1, j
        p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
        h(i,k) = h(i,k) - p
        h(i,k+1) = h(i,k+1) - p * q
        h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue
!
  260 continue
!
  go to 70
!    ........ one root found.
  270 wr(en) = x + t
  wi(en) = 0.0
  en = na
  go to 60
!    ........ two roots found.
  280 p = (y - x) / 2.0
  q = p * p + w
  zz = sqrt(abs(q))
  x = x + t
  if (q < 0.0) go to 320
!    ........ real pair.
  zz = p + sign(zz,p)
  wr(na) = x + zz
  wr(en) = wr(na)
  if (zz /= 0.0) wr(en) = x - w / zz
  wi(na) = 0.0
  wi(en) = 0.0
  go to 330
!    ........ complex pair.
  320 wr(na) = x + p
  wr(en) = x + p
  wi(na) = zz
  wi(en) = -zz
  330 en = enm2
  go to 60
!    ........ set error -- all eigenvalues have not
!                converged after 30*n iterations.
 1000 ierr = en
 1001 return
end
subroutine ibbs (ldd,ldt,n,kblszz,nsize,lbhb,iblock,d,t,jt,x,ivers,wksp)
!
!*******************************************************************************
!
!! IBBS does an incomplete block backward pass.
!
!
!     symmetric diagonal data structure, natural ordering.
!     block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         kblsz    block size
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         lbhb     number of blocks per block row
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         x        input/output vector of length n
!         ivers    key for version of factorization
!                   = 1   version 1
!                   = 2   version 2
!         wksp     real workspace vector
!
  integer   jt(1), iblock(3,1)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   vers2
!
  kblsz = kblszz
  l = n/kblsz
  nt = iblock(3,1) - 1
  vers2 = ivers == 2
  do 40 k = l,1,-1
     ist = (k - 1)*kblsz + 1
     ied = k*kblsz
     if (k == l) go to 15
     jjlim = min (lbhb,l-k+2)
     do 10 jj = 3,jjlim
        jblk = iblock(1,jj)
        jst = iblock(2,jj)
        mjj = iblock(3,jj)
        inc = jblk*kblsz
        istf = ist + inc
        if (istf > n) go to 10
        call vsubd (ldt,1,kblsz,kblsz,mjj,t(ist,jst),jt(jst), &
                       x(ist),x(istf),inc)
 10      continue
 15      if (nt >= 1) go to 25
     do 20 i = ist,ied
 20      x(i) = d(i,1)*x(i)
     go to 40
 25      if (vers2) go to 30
     call bdsol (ldd,kblsz,nsize,nt,0,d(ist,1),x(ist),x(ist),0)
     go to 40
 30      call bmul (ldd,kblsz,nt,d(ist,1),d(ist,2),x(ist),wksp)
     do 35 i = ist,ied
 35      x(i) = wksp(i-ist+1)
 40   continue
  return
end
subroutine ibbsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBBSN does an incomplete block backward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   unif, vers2
!
  vers2 = ivers == 2
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do backward solution.
!
 10   lm1 = l - 1
  do 50 k = lm1,1,-1
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     do 22 i = 1,na
 22      wksp(i) = 0.0
     do 25 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol <= k) go to 25
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        if (istb > n) go to 25
        call vaddd (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb),wksp, &
          x(istb),inc)
 25      continue
     if (ndt + ndb >= 1) go to 35
     do 30 i = ist,ied
 30      x(i) = x(i) - d(i,1)*wksp(i-ist+1)
     go to 50
 35      if (vers2) go to 40
     call bdsol (ldd,na,nsize,ndt,ndb,d(ist,1),wksp,wksp,1)
     do 37 i = ist,ied
 37      x(i) = x(i) - wksp(i-ist+1)
     go to 50
 40      nap1 = na + 1
     call bmuln (ldd,na,ndt,ndb,d(ist,1),d(ist,2),d(ist,ndt+2),wksp,wksp(nap1))
     do 45 i = ist,ied
 45      x(i) = x(i) - wksp(i-ist+nap1)
 50   continue
  return
end
subroutine ibbsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBBSNT does an incomplete block transpose backward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   unif, vers1
!
  vers1 = ivers == 1
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do backward solution.
!
 10   do 45 k = l,1,-1
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     if (ndt + ndb >= 1) go to 30
     do 25 i = ist,ied
 25      x(i) = d(i,1)*x(i)
     go to 35
 30      if (vers1) call bdsolt(ldd,na,nsize,ndt,ndb,d(ist,1),x(ist),x(ist))
     if (vers1) go to 35
     call bmulnt (ldd,na,ndt,ndb,d(ist,1),d(ist,2),d(ist,ndt+2),x(ist),wksp)

     do 32 i = ist,ied
 32      x(i) = wksp(i-ist+1)

 35      do 40 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol >= k) go to 40
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        if (istb < 1) go to 40
        call vsubdt (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb), &
                      x(istb),x(ist),inc)
 40      continue
 45   continue
  return
end
subroutine ibfcn1 (lddd,ldtt,n,jd,jt,d,t,ncol,nci,iblock,lbhb,iunif,ipropa, &
  ipt,omega,wksp,ier)
!
!*******************************************************************************
!
!! IBFCN1 does an incomplete block factorization.
!
!
!  The matrix is contained in d and t (version 1, unmodified).
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic (version 1) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer array of size ncolor by whatever
!                   giving the diagonal block diagonal numbers for
!                   each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jd(ncol,1), jt(ncol,1), nci(1), lbhb(1), iblock(3,ncol,2)
  dimension d(lddd,1), t(ldtt,1), wksp(1)
  logical   unif, propa
!
  ldd = lddd
  ldt = ldtt
  ncolor = ncol
  unif = iunif == 1
  propa = ipropa == 1
!
!  define various constants.
!
  if (unif) go to 15
  klim = ncolor
  go to 20
 15   kblsz = nci(1)
  na = kblsz
  nb = kblsz
  nc = kblsz
  ii = 1
  kk = 1
  jlim = lbhb(1)
  llim = jlim
  klim = n/kblsz
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  ma = ndt + ndb + 1
!
!  start factorization.
!
 20   do 95 k = 1,klim
     if (unif) go to 25
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     ma = ndt + ndb + 1
     go to 30
 25      ist = (k - 1)*kblsz + 1
 30      call bdfac (ldd,na,na,ndt,ndb,d(ist,1),1)
     call mcopy (ldd,na,na,ma,d(ist,1),wksp)
     call bdinv (na,na,na,ndt,ndb,wksp,1)
     if (k == klim .or. jlim <= 2) go to 95
     do 90 i = k+1,klim
        if (unif) go to 35
        ii = i
        llim = lbhb(i)
 35         if (llim <= 2) go to 90
        do 40 l = 3,llim
           jcol = i + iblock(1,ii,l)
           if (jcol == k) go to 45
 40         continue
        go to 90
 45         mc = iblock(3,ii,l)
        if (unif) go to 50
        nc = ipt(i+1) - ipt(i)
        incc = ipt(k) - ipt(i)
        go to 55
 50         incc = (k - i)*kblsz
 55         istc = ist - incc
        jstc = iblock(2,ii,l)
        do 85 j = 3,jlim
           jcol = k + iblock(1,kk,j)
           if (jcol <= k) go to 85
           jdiff = jcol - i
           if (jdiff /= 0 .and. propa) go to 85
           do 60 m = 1,llim
              if (iblock(1,ii,m) == jdiff) go to 65
 60            continue
           go to 85
 65            mb = iblock(3,kk,j)
           istb = ist
           jstb = iblock(2,kk,j)
           if (unif) go to 70
           nb = ipt(jcol+1) - ipt(jcol)
           incb = ipt(jcol) - ipt(k)
           go to 75
 70            incb = (jcol - k)*kblsz
 75            incd = incc + incb
           istd = istc
           jstd = iblock(2,ii,m)
           md = iblock(3,ii,m)
           if (m == 1) go to 80
           call t1prod (na,ldt,ldt,ldt,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jt(ii,jstd),wksp,t(istb,jstb), &
                           t(istc,jstc),t(istd,jstd))
           go to 85
 80            md = md + iblock(3,ii,2)
           call t1prod (na,ldt,ldt,ldd,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jd(ii,jstd),wksp,t(istb,jstb), &
                           t(istc,jstc),d(istd,jstd))
 85         continue
 90      continue
 95   continue
  return
end
subroutine ibfcn2 (lddd,ldtt,n,jd,jt,d,t,ncol,nci,iblock,lbhb,iunif, &
  ipropa,ipt,omega,wksp,ier)
!
!*******************************************************************************
!
!! IBFCN2 does an incomplete block factorization.
!
!  the matrix is contained in d and t (version 2, unmodified).
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic (version 2) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer array of size ncolor by whatever
!                   giving the diagonal block diagonal numbers for
!                   each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0
!
!  
!
  integer   ipt(1), jd(ncol,1), jt(ncol,1), nci(1), lbhb(1), iblock(3,ncol,2)
  dimension d(lddd,1), t(ldtt,1), wksp(1)
  logical   unif, propa
!
  ldd = lddd
  ldt = ldtt
  ncolor = ncol
  unif = iunif == 1
  propa = ipropa == 1
!
!  define various constants.
!
  if (unif) go to 15
  klim = ncolor
  go to 20
 15   kblsz = nci(1)
  na = kblsz
  nb = kblsz
  nc = kblsz
  ii = 1
  kk = 1
  jlim = lbhb(1)
  llim = jlim
  klim = n/kblsz
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  ma = ndt + ndb + 1
!
!  start factorization.
!
 20   do 95 k = 1,klim
     if (unif) go to 25
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     ma = ndt + ndb + 1
     go to 30
 25      ist = (k - 1)*kblsz + 1
 30      call bdfac (ldd,na,na,ndt,ndb,d(ist,1),1)
     call bdinv (ldd,na,na,ndt,ndb,d(ist,1),1)
     if (k == klim .or. jlim <= 2) go to 95
     do 90 i = k+1,klim
        if (unif) go to 35
        ii = i
        llim = lbhb(i)
 35         if (llim <= 2) go to 90
        do 40 l = 3,llim
           jcol = i + iblock(1,ii,l)
           if (jcol == k) go to 45
 40         continue
        go to 90
 45         mc = iblock(3,ii,l)
        if (unif) go to 50
        nc = ipt(i+1) - ipt(i)
        incc = ipt(k) - ipt(i)
        go to 55
 50         incc = (k - i)*kblsz
 55         istc = ist - incc
        jstc = iblock(2,ii,l)
        do 85 j = 3,jlim
           jcol = k + iblock(1,kk,j)
           if (jcol <= k) go to 85
           jdiff = jcol - i
           if (jdiff /= 0 .and. propa) go to 85
           do 60 m = 1,llim
              if (iblock(1,ii,m) == jdiff) go to 65
 60            continue
           go to 85
 65            mb = iblock(3,kk,j)
           istb = ist
           jstb = iblock(2,kk,j)
           if (unif) go to 70
           nb = ipt(jcol+1) - ipt(jcol)
           incb = ipt(jcol) - ipt(k)
           go to 75
 70            incb = (jcol - k)*kblsz
 75            incd = incc + incb
           istd = istc
           jstd = iblock(2,ii,m)
           md = iblock(3,ii,m)
           if (m == 1) go to 80
           call t1prod (ldd,ldt,ldt,ldt,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jt(ii,jstd),d(ist,1),t(istb,jstb), &
                           t(istc,jstc),t(istd,jstd))
           go to 85
 80            md = md + iblock(3,ii,2)
           call t1prod (ldd,ldt,ldt,ldd,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jd(ii,jstd),d(ist,1),t(istb,jstb), &
                           t(istc,jstc),d(istd,jstd))
 85         continue
 90      continue
 95   continue
  return
end
subroutine ibfcn3 (lddd,ldtt,n,jd,jt,d,t,ncol,nci,iblock,lbhb,iunif, &
  ipropa,ipt,omega,wksp,ier)
!
!*******************************************************************************
!
!! IBFCN3 does an incomplete block factorization.
!
!
!  The matrix is contained in d and t (version 1, modified).
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic (version 1) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer array of size ncolor by whatever
!                   giving the diagonal block diagonal numbers for
!                   each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0
!         omega    relaxation factor between 0 and 1.
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jd(ncol,1), jt(ncol,1), nci(1), lbhb(1), iblock(3,ncol,2)
  dimension d(lddd,1), t(ldtt,1), wksp(1)
  logical   unif, propa
!
  ldd = lddd
  ldt = ldtt
  ncolor = ncol
  unif = iunif == 1
  propa = ipropa == 1
!
!  define various constants.
!
  if (unif) go to 15
  klim = ncolor
  go to 20
 15   kblsz = nci(1)
  na = kblsz
  nb = kblsz
  nc = kblsz
  ii = 1
  kk = 1
  jlim = lbhb(1)
  llim = jlim
  klim = n/kblsz
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  ma = ndt + ndb + 1
!
!  start factorization.
!
 20   do 100 k = 1,klim
     if (unif) go to 25
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     ma = ndt + ndb + 1
     go to 30
 25      ist = (k - 1)*kblsz + 1
 30      call bdfac (ldd,na,na,ndt,ndb,d(ist,1),1)
     call mcopy (ldd,na,na,ma,d(ist,1),wksp)
     call bdinv (na,na,na,ndt,ndb,wksp,1)
     ip1 = na*ma + 1
     ip2 = ip1 + na - 1
     if (k == klim .or. jlim <= 2) go to 100
     do 95 i = k+1,klim
        if (unif) go to 35
        ii = i
        llim = lbhb(i)
 35         if (llim <= 2) go to 95
        do 40 l = 3,llim
           jcol = i + iblock(1,ii,l)
           if (jcol == k) go to 45
 40         continue
        go to 95
 45         mc = iblock(3,ii,l)
        if (unif) go to 50
        nc = ipt(i+1) - ipt(i)
        incc = ipt(k) - ipt(i)
        go to 55
 50         incc = (k - i)*kblsz
 55         istc = ist - incc
        jstc = iblock(2,ii,l)
        do 90 j = 3,jlim
           jcol = k + iblock(1,kk,j)
           if (jcol <= k) go to 90
           mb = iblock(3,kk,j)
           istb = ist
           jstb = iblock(2,kk,j)
           if (unif) go to 60
           nb = ipt(jcol+1) - ipt(jcol)
           incb = ipt(jcol) - ipt(k)
           go to 65
 60            incb = (jcol - k)*kblsz
 65            incd = incc + incb
           istd = istc
           jdiff = jcol - i
           if (jdiff /= 0 .and. propa) go to 85
           do 70 m = 1,llim
              if (iblock(1,ii,m) == jdiff) go to 75
 70            continue
           go to 85
 75            jstd = iblock(2,ii,m)
           md = iblock(3,ii,m)
           if (m == 1) go to 80
           call t1prod (na,ldt,ldt,ldt,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jt(ii,jstd),wksp,t(istb,jstb), &
                           t(istc,jstc),t(istd,jstd))
           call tsumn &
                    (na,nc,nb,na,ldt,ldt,ncolor,ma,mb,mc,md,incb, &
                     incc,incd,jd(kk,1),jt(kk,jstb),jt(ii,jstc), &
                     jt(ii,jstd),wksp,t(istb,jstb),t(istc,jstc), &
                     d(istd,1),omega)
           go to 85
 80            md = md + iblock(3,ii,2)
           call t1prod (na,ldt,ldt,ldd,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jd(ii,jstd),wksp,t(istb,jstb), &
                           t(istc,jstc),d(istd,jstd))
           call tsumn &
                    (na,nc,nb,na,ldt,ldt,ncolor,ma,mb,mc,md,incb, &
                     incc,incd,jd(kk,1),jt(kk,jstb),jt(ii,jstc), &
                     jd(ii,jstd),wksp,t(istb,jstb),t(istc,jstc), &
                     d(istd,1),omega)
 85            call rowsum (ldt,na,mb,t(istb,jstb),wksp(ip1),1)
           do 87 iii = ip1,ip2
 87            wksp(iii) = omega*wksp(iii)
           call bdsol (ldd,na,na,ndt,ndb,d(ist,1),wksp(ip1),wksp(ip1),1)
           call vsubd (ldt,ncolor,nc,na,mc,t(istc,jstc), &
                           jt(ii,jstc),d(istd,1),wksp(ip1),incc)
 90         continue
 95      continue
 100  continue
  return
end
subroutine ibfcn4 (lddd,ldtt,n,jd,jt,d,t,ncol,nci,iblock,lbhb,iunif, &
  ipropa,ipt,omega,wksp,ier)
!
!*******************************************************************************
!
!! IBFCN4 does an incomplete block factorization.
!
!
!   The matrix is contained in d and t (version 2, modified).
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic (version 2) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer array of size ncolor by whatever
!                   giving the diagonal block diagonal numbers for
!                   each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0
!         omega    relaxation factor between 0 and 1.
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jd(ncol,1), jt(ncol,1), nci(1), lbhb(1), iblock(3,ncol,2)
  dimension d(lddd,2), t(ldtt,1), wksp(1)
  logical   unif, propa
!
  ldd = lddd
  ldt = ldtt
  ncolor = ncol
  unif = iunif == 1
  propa = ipropa == 1
!
!  define various constants.
!
  ip1 = n + 1
  if (unif) go to 15
  klim = ncolor
  do 13 k = 1,ncolor
     ist = ipt(k) + 1
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     ma = ndt + ndb + 1
     call rowsum (ldd,na,ma,d(ist,1),wksp(ist),1)
 13   continue
  go to 20
 15   kblsz = nci(1)
  na = kblsz
  nb = kblsz
  nc = kblsz
  ii = 1
  kk = 1
  jlim = lbhb(1)
  llim = jlim
  klim = n/kblsz
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  ma = ndt + ndb + 1
  call rowsum (ldd,n,ma,d,wksp,1)
!
!  start factorization.
!
 20   do 100 k = 1,klim
     if (unif) go to 25
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     ma = ndt + ndb + 1
     go to 30
 25      ist = (k - 1)*kblsz + 1
 30      isu = ist + na - 1
     call bdfac (ldd,na,na,ndt,ndb,d(ist,1),1)
     call bdinv (ldd,na,na,ndt,ndb,d(ist,1),1)
     call bmuln (ldd,na,ndt,ndb,d(ist,1),d(ist,2),d(ist,ndt+2), &
                    wksp(ist),wksp(ip1))
     do 31 iii = ist,isu
        if (wksp(iii) /= 0.0) go to 31
        ier = -12
        call ershow (ier,'ibfcn4')
        return
 31      continue
     do 33 iii = ist,isu
 33      d(iii,1) = d(iii,1) + omega*(1.0 - wksp(iii-ist+ip1))/wksp(iii)
     ip2 = ip1 + na
     if (k == klim .or. jlim <= 2) go to 100
     do 95 i = k+1,klim
        if (unif) go to 35
        ii = i
        llim = lbhb(i)
 35         if (llim <= 2) go to 95
        do 40 l = 3,llim
           jcol = i + iblock(1,ii,l)
           if (jcol == k) go to 45
 40         continue
        go to 95
 45         mc = iblock(3,ii,l)
        if (unif) go to 50
        nc = ipt(i+1) - ipt(i)
        incc = ipt(k) - ipt(i)
        go to 55
 50         incc = (k - i)*kblsz
 55         istc = ist - incc
        jstc = iblock(2,ii,l)
        do 90 j = 3,jlim
           jcol = k + iblock(1,kk,j)
           if (jcol <= k) go to 90
           mb = iblock(3,kk,j)
           istb = ist
           jstb = iblock(2,kk,j)
           if (unif) go to 60
           nb = ipt(jcol+1) - ipt(jcol)
           incb = ipt(jcol) - ipt(k)
           go to 65
 60            incb = (jcol - k)*kblsz
 65            incd = incc + incb
           istd = istc
           jdiff = jcol - i
           if (jdiff /= 0 .and. propa) go to 85
           do 70 m = 1,llim
              if (iblock(1,ii,m) == jdiff) go to 75
 70            continue
           go to 85
 75            jstd = iblock(2,ii,m)
           md = iblock(3,ii,m)
           if (m == 1) go to 80
           call t1prod (ldd,ldt,ldt,ldt,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jt(ii,jstd),d(ist,1),t(istb,jstb), &
                           t(istc,jstc),t(istd,jstd))
           call tsumn &
                    (na,nc,nb,ldd,ldt,ldt,ncolor,ma,mb,mc,md,incb, &
                     incc,incd,jd(kk,1),jt(kk,jstb),jt(ii,jstc), &
                     jt(ii,jstd),d(ist,1),t(istb,jstb),t(istc,jstc), &
                     wksp(istd),1.0)
           go to 85
 80            md = md + iblock(3,ii,2)
           call t1prod (ldd,ldt,ldt,ldd,ncolor,na,nc,nb, &
                           ma,mb,mc,md,incb,incc,incd,jd(kk,1), &
                           jt(kk,jstb),jt(ii,jstc), &
                           jd(ii,jstd),d(ist,1),t(istb,jstb), &
                           t(istc,jstc),d(istd,jstd))
 85            call rowsum (ldt,na,mb,t(istb,jstb),wksp(ip1),1)
           call bmuln (ldd,na,ndt,ndb,d(ist,1),d(ist,2), &
                          d(ist,ndt+2),wksp(ip1),wksp(ip2))
           call vsubd (ldt,ncolor,nc,na,mc,t(istc,jstc), &
                           jt(ii,jstc),wksp(istd),wksp(ip2),incc)
 90         continue
 95      continue
 100  continue
  return
end
subroutine ibfcs1 (lddd,ldtt,nn,jd,jt,d,t,kblszz, iblock,lbhb,ipropa,omega, &
  wksp,ier)
!
!*******************************************************************************
!
!! IBFCS1 does an incomplete block factorization.
!
!
!  The matrix is contained in d and t (version 1, unmodified).
!     symmetric diagonal data structure, natural ordering.
!     block ic (version 1) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer vector giving the diagonal numbers
!                   for the diagonal block
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         kblsz    block size
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         lbhb     number of blocks per block row
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!         wksp     real workspace vector
!
!  
!
  integer   jd(1), jt(1), iblock(3,3)
  dimension d(lddd,1), t(ldtt,1), wksp(1)
  logical   propa
!
  n = nn
  ldd = lddd
  ldt = ldtt
  na = kblszz
  propa = ipropa == 1
  klim = n/na
  ma = iblock(3,1)
  ndt = ma - 1
!
!  block tridiagonal case.
!
  if (lbhb > 3) go to 25
  jblkb = iblock(1,3)
  mb = iblock(3,3)
  incb = jblkb*na
  do 20 k = 1,klim
     ist = (k - 1)*na + 1
     istd = ist + incb
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     if (istd > n) go to 20
     call mcopy (ldd,na,na,ma,d(ist,1),wksp)
     call bdinv (na,na,na,ndt,0,wksp,0)
     call t2prod (na,na,ldt,ldt,ldd,ma,mb,mb,ma,incb,incb,0,jd,jt,jt,jd, &
       wksp,t(ist,1),t(ist,1),d(istd,1))
 20   continue
  return
!
!  general block structure.
!
 25   do 50 k = 1,klim
     ist = (k - 1)*na + 1
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     if (k == klim) go to 50
     call mcopy (ldd,na,na,ma,d(ist,1),wksp)
     call bdinv (na,na,na,ndt,0,wksp,0)
     jjlim = min(lbhb,klim-k+2)
     do 45 jjc = 3,jjlim
        jblkc = iblock(1,jjc)
        jstc = iblock(2,jjc)
        mc = iblock(3,jjc)
        incc = jblkc*na
        istd = ist + incc
        if (istd > n) go to 45
        do 40 jjb = 3,jjlim
           jblkb = iblock(1,jjb)
           jstb = iblock(2,jjb)
           mb = iblock(3,jjb)
           incb = jblkb*na
           jdiff = jblkb - jblkc
           if (jdiff < 0) go to 40
           if (jdiff /= 0 .and. propa) go to 40
           do 30 jjd = 1,jjlim
              if (jdiff == iblock(1,jjd)) go to 35
 30            continue
           go to 40
 35            jblkd = iblock(1,jjd)
           jstd = iblock(2,jjd)
           md = iblock(3,jjd)
           incd = jblkd*na
           if (jjd /= 1) call t2prod(na,na,ldt,ldt,ldt,ma,mb,mc,md,incb, &
             incc,incd,jd,jt(jstb),jt(jstc), &
             jt(jstd),wksp,t(ist,jstb),t(ist,jstc),t(istd,jstd))
           if (jjd == 1) call t2prod(na,na,ldt,ldt,ldd,ma,mb,mc,md,incb, &
             incc,incd,jd,jt(jstb),jt(jstc), &
             jd,wksp,t(ist,jstb),t(ist,jstc),d(istd,1))
 40         continue
 45      continue
 50   continue
  return
end
subroutine ibfcs2 (lddd,ldtt,nn,jd,jt,d,t,kblszz,iblock,lbhb,ipropa, &
  omega,wksp,ier)
!
!*******************************************************************************
!
!! IBFCS2 does an incomplete block factorization.
!
!
!  The matrix is contained in d and t (version 2, unmodified).
!     symmetric diagonal data structure, natural ordering.
!     block ic (version 2) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer vector giving the diagonal numbers
!                   for the diagonal block
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         kblsz    block size
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         lbhb     number of blocks per block row
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!
!  
!
  integer   jd(1), jt(1), iblock(3,3)
  dimension d(lddd,1), t(ldtt,1), wksp(1)
  logical   propa
!
  n = nn
  ldd = lddd
  ldt = ldtt
  na = kblszz
  propa = ipropa == 1
  klim = n/na
  ma = iblock(3,1)
  ndt = ma - 1
!
!  block tridiagonal case.
!
  if (lbhb > 3) go to 25
  jblkb = iblock(1,3)
  mb = iblock(3,3)
  incb = jblkb*na
  do 20 k = 1,klim
     ist = (k - 1)*na + 1
     istd = ist + incb
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     call bdinv (ldd,na,na,ndt,0,d(ist,1),0)
     if (istd > n) go to 20
     call t2prod (na,ldd,ldt,ldt,ldd,ma,mb,mb,ma,incb,incb,0,jd,jt,jt, &
       jd,d(ist,1),t(ist,1),t(ist,1),d(istd,1))
 20   continue
  return
!
!  general block structure.
!
 25   do 50 k = 1,klim
     ist = (k - 1)*na + 1
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     call bdinv (ldd,na,na,ndt,0,d(ist,1),0)
     if (k == klim) go to 50
     jjlim = min(lbhb,klim-k+2)
     do 45 jjc = 3,jjlim
        jblkc = iblock(1,jjc)
        jstc = iblock(2,jjc)
        mc = iblock(3,jjc)
        incc = jblkc*na
        istd = ist + incc
        if (istd > n) go to 45
        do 40 jjb = 3,jjlim
           jblkb = iblock(1,jjb)
           jstb = iblock(2,jjb)
           mb = iblock(3,jjb)
           incb = jblkb*na
           jdiff = jblkb - jblkc
           if (jdiff < 0) go to 40
           if (jdiff /= 0 .and. propa) go to 40
           do 30 jjd = 1,jjlim
              if (jdiff == iblock(1,jjd)) go to 35
 30            continue
           go to 40
 35            jblkd = iblock(1,jjd)
           jstd = iblock(2,jjd)
           md = iblock(3,jjd)
           incd = jblkd*na
           if (jjd /= 1) call t2prod(na,ldd,ldt,ldt,ldt,ma,mb,mc,md,incb, &
             incc,incd,jd,jt(jstb),jt(jstc), &
             jt(jstd),d(ist,1),t(ist,jstb),t(ist,jstc),t(istd,jstd))
           if (jjd == 1) call t2prod(na,ldd,ldt,ldt,ldd,ma,mb,mc,md,incb, &
             incc,incd,jd,jt(jstb),jt(jstc), &
             jd,d(ist,1),t(ist,jstb),t(ist,jstc),d(istd,1))
 40         continue
 45      continue
 50   continue
  return
end
subroutine ibfcs3 (lddd,ldtt,nn,jd,jt,d,t,kblszz,iblock,lbhb,ipropa, &
  omegaa,wksp,ier)
!
!*******************************************************************************
!
!! IBFCS3 does an incomplete block factorization.
!
!
!  The matrix is contained in d and t (version 1, modified).
!     symmetric diagonal data structure, natural ordering.
!     block ic (version 1) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer vector giving the diagonal numbers
!                   for the diagonal block
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         kblsz    block size
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         lbhb     number of blocks per block row
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!         omega    relaxation factor between 0. and 1.
!                   = 0   no modification
!                   = 1   full modification
!         wksp     real workspace vector
!
!  
!
  integer   jd(1), jt(1), iblock(3,3)
  dimension d(lddd,1), t(ldtt,1), wksp(1)
  logical   propa
!
  n = nn
  ldd = lddd
  ldt = ldtt
  na = kblszz
  omega = omegaa
  propa = ipropa == 1
  klim = n/na
  ma = iblock(3,1)
  ndt = ma - 1
!
!  block tridiagonal case.
!
  if (lbhb > 3) go to 25
  ip1 = na*ma + 1
  ip2 = ip1 + na - 1
  jblkb = iblock(1,3)
  mb = iblock(3,3)
  incb = jblkb*na
  do 20 k = 1,klim
     ist = (k - 1)*na + 1
     istd = ist + incb
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     if (istd > n) go to 20
     call mcopy (ldd,na,na,ma,d(ist,1),wksp)
     call bdinv (na,na,na,ndt,0,wksp,0)
     call t2prod (na,na,ldt,ldt,ldd,ma,mb,mb,ma,incb,incb,0,jd,jt,jt,jd, &
       wksp,t(ist,1),t(ist,1),d(istd,1))
     call tsum (na,na,ldt,ldt,ma,mb,mb,ma,incb,incb,0,jd,jt,jt,jd,wksp, &
       t(ist,1),t(ist,1),d(istd,1),d(istd,1),wksp(ip1),1,omega)
     call rowsum (ldt,na,mb,t(ist,1),wksp(ip1),1)
     do 15 iii = ip1,ip2
 15      wksp(iii) = omega*wksp(iii)
     call bdsol (ldd,na,na,ndt,0,d(ist,1),wksp(ip1),wksp(ip1),0)
     call vsubdt (ldt,1,na,na,mb,t(ist,1),jt,d(istd,1),wksp(ip1),incb)
 20   continue
  return
!
!  general block structure.
!
 25   ip1 = na*ma + 1
  ip2 = ip1 + na - 1
  do 60 k = 1,klim
     ist = (k - 1)*na + 1
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     if (k == klim) go to 60
     call mcopy (ldd,na,na,ma,d(ist,1),wksp)
     call bdinv (na,na,na,ndt,0,wksp,0)
     jjlim = min(lbhb,klim-k+2)
     do 55 jjc = 3,jjlim
        jblkc = iblock(1,jjc)
        jstc = iblock(2,jjc)
        mc = iblock(3,jjc)
        incc = jblkc*na
        istd = ist + incc
        if (istd > n) go to 55
        do 50 jjb = 3,jjlim
           jblkb = iblock(1,jjb)
           jstb = iblock(2,jjb)
           mb = iblock(3,jjb)
           incb = jblkb*na
           istdd = ist + incb
           if (istdd > n) go to 50
           jdiff = jblkb - jblkc
           if (jdiff < 0) go to 50
           if (jdiff /= 0 .and. propa) go to 40
           do 30 jjd = 1,jjlim
              if (jdiff == iblock(1,jjd)) go to 35
 30            continue
           go to 40
 35            jblkd = iblock(1,jjd)
           jstd = iblock(2,jjd)
           md = iblock(3,jjd)
           incd = jblkd*na
           if (jjd /= 1) call t2prod(na,na,ldt,ldt,ldt,ma,mb,mc,md,incb, &
                     incc,incd,jd,jt(jstb),jt(jstc), &
                     jt(jstd),wksp,t(ist,jstb),t(ist,jstc), &
                     t(istd,jstd))
           if (jjd == 1) call t2prod(na,na,ldt,ldt,ldd,ma,mb,mc,md,incb, &
                     incc,incd,jd,jt(jstb),jt(jstc), &
                     jd,wksp,t(ist,jstb),t(ist,jstc), &
                     d(istd,1))
           if (jjd /= 1) call tsum(na,na,ldt,ldt,ma,mb,mc,md,incb, &
                     incc,incd,jd,jt(jstb),jt(jstc), &
                     jt(jstd),wksp,t(ist,jstb),t(ist,jstc), &
                     d(istd,1),d(istdd,1),wksp(ip1),0,omega)
           if (jjd == 1) call tsum(na,na,ldt,ldt,ma,mb,mc,md,incb, &
                     incc,incd,jd,jt(jstb),jt(jstc), &
                     jd,wksp,t(ist,jstb),t(ist,jstc), &
                     d(istd,1),d(istdd,1),wksp(ip1),1,omega)
!
 40            call rowsum (ldt,na,mb,t(ist,jstb),wksp(ip1),1)
           do 42 iii = ip1,ip2
 42            wksp(iii) = omega*wksp(iii)
           call bdsol (ldd,na,na,ndt,0,d(ist,1),wksp(ip1),wksp(ip1),0)
           call vsubdt (ldt,1,na,na,mc,t(ist,jstc),jt(jstc), &
                         d(istd,1),wksp(ip1),incc)
           if (jdiff == 0) go to 50
           call rowsum (ldt,na,mc,t(ist,jstc),wksp(ip1),1)
           do 45 iii = ip1,ip2
 45            wksp(iii) = omega*wksp(iii)
           call bdsol (ldd,na,na,ndt,0,d(ist,1),wksp(ip1),wksp(ip1),0)
           call vsubdt (ldt,1,na,na,mb,t(ist,jstb),jt(jstb), &
                       d(istdd,1),wksp(ip1),incb)
 50         continue
 55      continue
 60   continue
  return
end
subroutine ibfcs4 (lddd,ldtt,nn,jd,jt,d,t,kblszz,iblock,lbhb,ipropa, &
  omegaa,wksp,ier)
!
!*******************************************************************************
!
!! IBFCS4 does an incomplete block factorization.
!
!
!  The matrix is contained in d and t (version 2, modified).
!     symmetric diagonal data structure, natural ordering.
!     block ic (version 2) preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         jd       integer vector giving the diagonal numbers
!                   for the diagonal block
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         kblsz    block size
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         lbhb     number of blocks per block row
!         ipropa   property a switch
!                   = 0   matrix does not have block property a
!                   = 1   matrix has block property a
!         omega    relaxation factor between 0. and 1.
!                   = 0   no modification
!                   = 1   full modification
!         wksp     real workspace vector
!
!  
!
  integer   jd(1), jt(1), iblock(3,3)
  dimension d(lddd,2), t(ldtt,1), wksp(1)
  logical   propa
!
  n = nn
  ldd = lddd
  ldt = ldtt
  na = kblszz
  omega = omegaa
  propa = ipropa == 1
  klim = n/na
  ma = iblock(3,1)
  ndt = ma - 1
!
!  block tridiagonal case.
!
  if (lbhb > 3) go to 25
  ip1 = n + 1
  ip2 = ip1 + na
  jblkb = iblock(1,3)
  mb = iblock(3,3)
  incb = jblkb*na
  call rowsum (ldd,n,ma,d,wksp,0)
  do 20 k = 1,klim
     ist = (k - 1)*na + 1
     isu = k*na
     istd = ist + incb
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     call bdinv (ldd,na,na,ndt,0,d(ist,1),0)
     call bmul (ldd,na,ndt,d(ist,1),d(ist,2),wksp(ist),wksp(ip1))
     do 10 ii = ist,isu
        if (wksp(ii) /= 0.0) go to 10
        ier = -12
        call ershow (ier,'ibfcs4')
        return
 10      continue
     do 15 ii = ist,isu
 15      d(ii,1) = d(ii,1) + omega*(1.0 - wksp(ii-ist+ip1))/wksp(ii)
     if (istd > n) go to 20
     call t2prod (na,ldd,ldt,ldt,ldd,ma,mb,mb,ma,incb,incb,0,jd,jt,jt,jd, &
       d(ist,1),t(ist,1),t(ist,1),d(istd,1))
     call rowsum (ldt,na,mb,t(ist,1),wksp(ip1),1)
     call bmul (ldd,na,ndt,d(ist,1),d(ist,2),wksp(ip1),wksp(ip2))
     call vsubdt (ldt,1,na,na,mb,t(ist,1),jt,wksp(istd),wksp(ip2),incb)
 20   continue
  return
!
!  general block structure.
!
 25   ip1 = n + 1
  ip2 = ip1 + na
  call rowsum (ldd,n,ma,d,wksp,0)
  do 60 k = 1,klim
     ist = (k - 1)*na + 1
     isu = k*na
     call bdfac (ldd,na,na,ndt,0,d(ist,1),0)
     call bdinv (ldd,na,na,ndt,0,d(ist,1),0)
     call bmul (ldd,na,ndt,d(ist,1),d(ist,2),wksp(ist),wksp(ip1))
     do 26 ii = ist,isu
        if (wksp(ii) /= 0.0) go to 26
        ier = -12
        call ershow (ier,'ibfcs4')
        return
 26      continue
     do 27 ii = ist,isu
 27      d(ii,1) = d(ii,1) + omega*(1.0 - wksp(ii-ist+ip1))/ wksp(ii)
     if (k == klim) go to 60
     jjlim = min(lbhb,klim-k+2)
     do 55 jjc = 3,jjlim
        jblkc = iblock(1,jjc)
        jstc = iblock(2,jjc)
        mc = iblock(3,jjc)
        incc = jblkc*na
        istd = ist + incc
        if (istd > n) go to 55
        do 50 jjb = 3,jjlim
           jblkb = iblock(1,jjb)
           jstb = iblock(2,jjb)
           mb = iblock(3,jjb)
           incb = jblkb*na
           istdd = ist + incb
           if (istdd > n) go to 50
           jdiff = jblkb - jblkc
           if (jdiff < 0) go to 50
           if (jdiff /= 0 .and. propa) go to 40
           do 30 jjd = 1,jjlim
              if (jdiff == iblock(1,jjd)) go to 35
 30            continue
           go to 40
 35            jblkd = iblock(1,jjd)
           jstd = iblock(2,jjd)
           md = iblock(3,jjd)
           incd = jblkd*na

           if (jjd /= 1) call t2prod(na,ldd,ldt,ldt,ldt,ma,mb,mc,md,incb, &
                     incc,incd,jd,jt(jstb),jt(jstc), &
                     jt(jstd),d(ist,1),t(ist,jstb),t(ist,jstc), &
                     t(istd,jstd))

           if (jjd == 1) call t2prod(na,ldd,ldt,ldt,ldd,ma,mb,mc,md,incb, &
                     incc,incd,jd,jt(jstb),jt(jstc), &
                     jd,d(ist,1),t(ist,jstb),t(ist,jstc), &
                     d(istd,1))

           if (jjd /= 1) call tsum(na,ldd,ldt,ldt,ma,mb,mc,md,incb, &
                     incc,incd,jd,jt(jstb),jt(jstc), &
                     jt(jstd),d(ist,1),t(ist,jstb),t(ist,jstc), &
                     wksp(istd),wksp(istdd),wksp(ip1),0,1.0)
!
 40            call rowsum (ldt,na,mb,t(ist,jstb),wksp(ip1),1)
           call bmul (ldd,na,ndt,d(ist,1),d(ist,2),wksp(ip1),wksp(ip2))
           call vsubdt (ldt,1,na,na,mc,t(ist,jstc),jt(jstc), &
                           wksp(istd),wksp(ip2),incc)
           if (jdiff == 0) go to 50
           call rowsum (ldt,na,mc,t(ist,jstc),wksp(ip1),1)
           call bmul (ldd,na,ndt,d(ist,1),d(ist,2),wksp(ip1),wksp(ip2))
           call vsubdt (ldt,1,na,na,mb,t(ist,jstb),jt(jstb), &
                           wksp(istdd),wksp(ip2),incb)
 50         continue
 55      continue
 60   continue
  return
end
subroutine ibfs (ldd,ldt,n,kblszz,nsize,lbhb,iblock,d,t,jt,x,ivers,wksp)
!
!*******************************************************************************
!
!! IBFS does an incomplete block forward pass.
!
!
!     symmetric diagonal data structure, natural ordering.
!     block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         kblsz    block size
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         lbhb     number of blocks per block row
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         x        input/output vector of length n
!         ivers    key for version of factorization
!                   = 1   version 1
!                   = 2   version 2
!         wksp     real workspace vector
!
!  
!
  integer   jt(1), iblock(3,1)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   vers1, vers2
!
  kblsz = kblszz
  l = n/kblsz
  lm1 = l - 1
  nt = iblock(3,1) - 1
  vers1 = ivers == 1
  vers2 = ivers == 2
  do 30 k = 1,lm1
     ist = (k - 1)*kblsz + 1
     ied = k*kblsz
     if (nt >= 1) go to 15
     do 10 i = ist,ied
 10      wksp(i-ist+1) = d(i,1)*x(i)
     go to 20
 15      if (vers1) call bdsol (ldd,kblsz,nsize,nt,0,d(ist,1),x(ist),wksp,0)
     if (vers2) call bmul (ldd,kblsz,nt,d(ist,1),d(ist,2),x(ist),wksp)
 20      jjlim = min (lbhb,l-k+2)
     do 25 jj = 3,jjlim
        jblk = iblock(1,jj)
        jst = iblock(2,jj)
        mjj = iblock(3,jj)
        inc = jblk*kblsz
        istf = ist + inc
        if (istf > n) go to 25
        call vsubdt (ldt,1,kblsz,kblsz,mjj,t(ist,jst),jt(jst),x(istf),wksp,inc)
 25      continue
 30   continue
  return
end
subroutine ibfsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBFSN does an incomplete block forward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1),iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   unif, vers2
!
  vers2 = ivers == 2
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do forward solution.
!
 10   do 50 k = 1,l
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     do 25 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol >= k) go to 25
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        if (istb < 1) go to 25
        call vsubd (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb),x(ist), &
          x(istb),inc)
 25      continue
     if (ndt + ndb >= 1) go to 35
     do 30 i = ist,ied
 30      x(i) = d(i,1)*x(i)
     go to 50
 35      if (vers2) go to 40
     call bdsol (ldd,na,nsize,ndt,ndb,d(ist,1),x(ist),x(ist),1)
     go to 50
 40      call bmuln (ldd,na,ndt,ndb,d(ist,1),d(ist,2),d(ist,ndt+2),x(ist),wksp)
     do 45 i = ist,ied
 45      x(i) = wksp(i-ist+1)
 50   continue
  return
end
subroutine ibfsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBFSNT does an incomplete block transpose forward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1),iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   unif, vers1, vers2
!
  vers1 = ivers == 1
  vers2 = ivers == 2
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do forward solution.
!
 10   lm1 = l - 1
  do 45 k = 1,lm1
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     if (ndt + ndb >= 1) go to 30
     do 25 i = ist,ied
 25      wksp(i-ist+1) = d(i,1)*x(i)
     go to 35
 30      if (vers1) call bdsolt(ldd,na,nsize,ndt,ndb,d(ist,1),x(ist),wksp)
     if (vers2) call bmulnt(ldd,na,ndt,ndb,d(ist,1),d(ist,2),d(ist,ndt+2), &
                   x(ist),wksp)
 35      do 40 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol <= k) go to 40
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        if (istb > n) go to 40
        call vsubdt (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb), &
                      x(istb),wksp,inc)
 40      continue
 45   continue
  return
end
subroutine ibsl (ldd,ldt,n,kblsz,nsize,lbhb,iblock,d,t,jt,y,x,ivers,wksp)
!
!*******************************************************************************
!
!! IBSL does an incomplete block solution.
!
!
!     symmetric diagonal data structure, natural ordering.
!     block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         kblsz    block size
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         lbhb     number of blocks per block row
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         y        input vector for the right-hand-side
!         x        output vector for the solution to q*x = y
!         ivers    key for version of factorization
!                   = 1   version 1
!                   = 2   version 2
!         wksp     real workspace vector
!
!  
!
  integer   jt(1), iblock(3,1)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call ibfs (ldd,ldt,n,kblsz,nsize,lbhb,iblock,d,t,jt,x,ivers,wksp)
  call ibbs (ldd,ldt,n,kblsz,nsize,lbhb,iblock,d,t,jt,x,ivers,wksp)
  return
end
subroutine ibsln (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBSLN does an incomplete block solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1),iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call ibfsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  call ibbsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  return
end
subroutine ibsln1 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBSLN1 does an incomplete block forward pass.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call ibfsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  return
end
subroutine ibsln2 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBSLN2 does an incomplete block backward pass.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call ibbsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  return
end
subroutine ibsln3 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBSLN3 does an incomplete block transpose back solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call ibbsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  return
end
subroutine ibsln4 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBSLN4 does an incomplete block transpose forward pass.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call ibfsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  return
end
subroutine ibslnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  ivers,iunif,wksp)
!
!*******************************************************************************
!
!! IBSLNT does an incomplete block transpose solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ic preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         ivers    key for version number
!                   = 1  version 1
!                   = 2  version 2
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call ibfsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  call ibbsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,ivers, &
    iunif,wksp)
  return
end
subroutine ic1 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp, iparm,rparm,ier)
!
!*******************************************************************************
!
!! IC1 drives the IC preconditioner.
!
!
  external accel, suba8, suba9, subq86, subq87, subq88
  external subq89, subq90, subq91, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
  n = nn
  if (ifact == 0 .and. lvfill > 0) go to 20
  call move1 (ndim,mdim,n,maxnz,jcoef,coef,maxt,maxb,ier)
  if (ier < 0) then
     call ershow (ier,'ic1')
     return
  end if
 20   t1 = timer (dummy)
  if (ifact == 1) call pfact1 (coef,jcoef,wksp,iwksp,n,1,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba8,suba9,subq86,subq87,subq88,subq89,subq90,subq91, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine ic2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! IC2 drives the symmetric IC preconditioner.
!
!
  external accel, suba1, subq13, subq14, subq15, subq16, subq17, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
!
  t1 = timer (dummy)
  if (ifact == 1) call pfact2 (coef,jcoef,wksp,iwksp,n,1,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  leniw = max (maxnz,nfacti)
  iwkpt1 = iipnt
  iipnt = iipnt + leniw
  call split (accel,suba1,suba1,subq13,subq13,subq14,subq15,subq16,subq17, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  iipnt = iipnt - leniw
  return
end
subroutine ic3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! IC3 drives the nonsymmetric IC preconditioner.
!
  external accel, suba4, suba5, subq48, subq49, subq50
  external subq51, subq52, subq53, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
!
  n = nn
  call needw ('ic3',1,iipnt,maxnz,ier)
  if (ier < 0) return
  call needw ('ic3',0,irpnt,n,ier)
  if (ier < 0) return
  if (ifact == 0 .and. lvfill > 0) go to 20
  call move2 (ndim,n,maxnz,jcoef,coef,wksp(irpnt),iwksp(iipnt),maxt,maxb)
 20   t1 = timer (dummy)
  if (ifact == 1) call pfact3 (coef,jcoef,wksp,iwksp,n,1,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  leniw = max (maxnz,nfacti)
  iwkpt1 = iipnt
  iipnt = iipnt + leniw
  call split (accel,suba4,suba5,subq48,subq49,subq50,subq51,subq52,subq53, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  iipnt = iipnt - leniw
  return
end
subroutine ic6 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! IC6 drives the IC preconditioner.
!
!     (multi-color ordering)
!
  external accel, suba8, suba9, sub104, sub105, sub106, sub107, sub108
  external sub109, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
  n = nn
  t1 = timer (dummy)
  if (ifact == 1) call pfactc (coef,jcoef,wksp,iwksp,n,1,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba8,suba9,sub104,sub105,sub106,sub107,sub108,sub109, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine icbs (ndim,nn,maxtt,jt,d,t,ipropa,irwise,iwksp,x)
!
!*******************************************************************************
!
!! ICBS does an IC back solve (natural ordering, diagonal storage).
!
!
!        (i + t)*x = y    if not property a
!        (i + d*t)*x = y  if property a
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        x      on input, x contains y
!               on output, x is the solution to back-solve
!
!  
!
  dimension x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
!
  logical   propa
!
  n = nn
  maxt = maxtt
  nm1 = n - 1
  propa = ipropa == 1
  if (maxt < 1) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 70
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxt
 15   iwksp(i) = n - jt(i)
!
!  determine nc, imax.
!
 20   nc = 1
  do 25 i = 1,maxt
     nterm = iwksp(i) + 1
     if (nterm <= nc) go to 25
     nc = nterm
     imax = i
 25   continue
  if (nc <= 1) return
  ndel = jt(imax)
  iend = nc - 1
  if (ndel > 1) go to 50
!
!  special case for first super diagonal.
!
  nc1 = 1
  do 30 i = 1,maxt
     if (i == imax) go to 30
     if (iwksp(i) > nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imax) = nc1 - 1
  if (propa) go to 40
  do 35 k = iend,nc1,-1
 35   x(k) = x(k) - t(k,imax)*x(k+1)
  go to 20
 40   do 45 k = iend,nc1,-1
 45   x(k) = x(k) - d(k)*t(k,imax)*x(k+1)
  go to 20
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imax) = iwksp(imax) - ndel
  ibeg = max (iend - ndel,0) + 1
  if (propa) go to 60
!dir$ ivdep
  do 55 i = ibeg,iend
 55   x(i) = x(i) - t(i,imax)*x(i+ndel)
  go to 20
!dir$ ivdep
 60   do 65 i = ibeg,iend
 65   x(i) = x(i) - d(i)*t(i,imax)*x(i+ndel)
  go to 20
!
!  rowwise algorithm.
!
 70   do 85 i = nm1,1,-1
     do 75 j = 1,maxt
 75      iwksp(j) = min (n,i+jt(j))
     sum = 0.0
     do 80 j = 1,maxt
 80      sum = sum + t(i,j)*x(iwksp(j))
     if (propa) sum = d(i)*sum
     x(i) = x(i) - sum
 85   continue
  return
end
subroutine icbscp (ndimr,ndimi,n,jc,d,c,ncolor,nc,nt,ipropa,wksp,x)
!
!*******************************************************************************
!
!! ICBSCP does a back IC solve.  (Purdue storage, multicolor)
!
!
!       (i + t)*x = y    if ipropa = 0
!       (d + t)*x = y    if ipropa = 1
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!          x      on input, x contains y
!                 on output, x is the solution to the back solve
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1)
  dimension d(1), c(ndimr,1), x(1), wksp(1)
  logical   propa
!
  propa = ipropa == 1
!
  ied = n
  do 25 icol = ncolor,1,-1
     npt = nc(icol)
     ist = ied - npt + 1
     j2 = nt(icol)
     call vsubp (ndimr,ndimi,npt,j2,c(ist,1),jc(ist,1),x(ist),x,wksp)
     if (.not. propa) go to 20
     do 15 i = ist,ied
 15      x(i) = x(i)*d(i)
 20      ied = ied - npt
 25   continue
  return
end
subroutine icbsct (ndimr,ndimi,n,jc,d,c,ncolor,nc,nt,nb,ipropa, wksp,x)
!
!*******************************************************************************
!
!! ICBSCT does a transpose back IC solve.  (Purdue storage, multicolor)
!
!
!     (i + (b**t))*x = y    if ipropa = 0
!     (d + (b**t))*x = y    if ipropa = 1
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length max(nc(i))
!          x      on input, x contains y
!                 on output, x is the solution to the back solve
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimr,1), x(1), wksp(1)
  logical   propa
!
  propa = ipropa == 1
!
  ied = n
  do 25 icol = ncolor,1,-1
     npt = nc(icol)
     ist = ied - npt + 1
     if (.not. propa) go to 20
     do 15 i = ist,ied
 15      x(i) = x(i)*d(i)
 20      j1 = nt(icol) + 1
     mj = nb(icol)
     call vsubpt (ndimr,ndimi,npt,mj,c(ist,j1),jc(ist,j1),x,x(ist), wksp)
     ied = ied - npt
 25   continue
  return
end
subroutine icbsp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
!
!*******************************************************************************
!
!! ICBSP does an IC back solve (natural ordering, Purdue storage).
!
!
!        (i + t)*x = y    if ipropa = 0
!        (d + t)*x = y    if ipropa = 1
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        x      on input, x contains y
!               on output, x is the solution to back-solve
!
!  
!
  dimension x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
  logical   propa
!
  propa = ipropa == 1
  if (maxt >= 1) go to 15
  if (.not. propa) return
  do 10 i = 1,n
 10   x(i) = x(i)*d(i)
  return
 15   do 25 i = n,1,-1
     sum = x(i)
     do 20 j = 1,maxt
        sum = sum - t(i,j)*x(jt(i,j))
 20      continue
     if (propa) sum = sum*d(i)
     x(i) = sum
 25   continue
  return
end
subroutine icbst (ndim,nn,maxbb,jb,d,b,ipropa,irwise,iwksp,x)
!
!*******************************************************************************
!
!! ICBST does an iC back solve (natural ordering, diagonal storage).
!
!
!        (i + (b**t))*x = y    if not property a
!        (i + d*(b**t))*x = y  if property a
!
!  Parameters:
!
!        ndim   row dimension of b array
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the sub-
!                diagonals of the factorization if not property a
!                or the sub-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxb
!        x      on input, x contains y
!               on output, x is the solution to back-solve
!
!  
!
  dimension x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
  logical   propa
!
  n = nn
  maxb = maxbb
  propa = ipropa == 1
  if (maxb < 1) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 70
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxb
 15   iwksp(i) = n + jb(i)
!
!  determine nc, imax.
!
 20   nc = 1
  do 25 i = 1,maxb
     nterm = iwksp(i) + 1
     if (nterm <= nc) go to 25
     nc = nterm
     imax = i
 25   continue
  if (nc <= 1) return
  ndel = -jb(imax)
  iend = nc - 1
  if (ndel > 1) go to 50
!
!  special case for first sub diagonal.
!
  nc1 = 1
  do 30 i = 1,maxb
     if (i == imax) go to 30
     if (iwksp(i) > nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imax) = nc1 - 1
  if (propa) go to 40
  do 35 k = iend,nc1,-1
 35   x(k) = x(k) - b(k+1,imax)*x(k+1)
  go to 20
 40   do 45 k = iend,nc1,-1
 45   x(k) = x(k) - d(k)*b(k+1,imax)*x(k+1)
  go to 20
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imax) = iwksp(imax) - ndel
  ibeg = max (iend - ndel,0) + 1
  if (propa) go to 60
!dir$ ivdep
  do 55 i = ibeg,iend
 55   x(i) = x(i) - b(i+ndel,imax)*x(i+ndel)
  go to 20
!dir$ ivdep
 60   do 65 i = ibeg,iend
 65   x(i) = x(i) - d(i)*b(i+ndel,imax)*x(i+ndel)
  go to 20
!
!  rowwise algorithm.
!
 70   if (propa) go to 90
  do 85 i = n,2,-1
     do 75 j = 1,maxb
 75      iwksp(j) = max (1,i+jb(j))
     term = x(i)
     do 80 j = 1,maxb
 80      x(iwksp(j)) = x(iwksp(j)) - b(i,j)*term
 85   continue
  return
 90   do 105 i = n,2,-1
     do 95 j = 1,maxb
 95      iwksp(j) = max (1,i+jb(j))
     term = x(i)
     do 100 j = 1,maxb
 100     x(iwksp(j)) = x(iwksp(j)) - d(iwksp(j))*b(i,j)*term
 105  continue
  return
end
subroutine icbstp (ndimr,ndimi,n,maxb,jb,d,b,ipropa,x)
!
!*******************************************************************************
!
!! ICBSTP does a transpose IC back solve (natural ordering, Purdue storage).
!
!
!        (i + (b**t))*x = y    if ipropa = 0
!        (d + (b**t))*x = y    if ipropa = 1
!
!  Parameters:
!
!        n      order of system
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the lower
!                triangle of the factorization if ipropa = 0
!                or the lower triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        x      on input, x contains y
!               on output, x is the solution to back-solve
!
!  
!
  dimension x(1), d(1), b(ndimr,1)
  integer   jb(ndimi,1)
  logical   propa
!
  propa = ipropa == 1
  if (maxb >= 1) go to 15
  if (.not. propa) return
  do 10 i = 1,n
 10   x(i) = x(i)*d(i)
  return
 15   do 25 i = n,1,-1
     if (propa) x(i) = x(i)*d(i)
     term = x(i)
     do 20 j = 1,maxb
        x(jb(i,j)) = x(jb(i,j)) - b(i,j)*term
 20      continue
 25   continue
  return
end
subroutine icf (ndim,nn,maxtt,jt,d,t,meth,ipropa,omega,wksp,iwksp,iflag)
!
!*******************************************************************************
!
!! ICF computes an incomplete factorization.  (symmetric diagonal storage)
!
!
!   The matrix is stored in d and t and the factorization replaces it.
!     
!
!  Parameters:
!
!          ndim   row dimension of t array
!          n      order of system (= nn)
!          maxt   number of columns in t array
!          jt     integer vector giving the diagonal indices of
!                  the corresponding columns in t
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          t      array of active size n by maxt giving the
!                  super-diagonals of the matrix
!          meth   point factorization wanted
!                  = 1   ic
!                  = 2   mic
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          omega  modification factor between 0.0 and 1.0
!                  = 0     no modification
!                  = 1     full modification
!          wksp   workspace vector of length n
!          iwksp  integer workspace of length maxt**2
!          iflag  indicator of factorization stability
!                    iflag = 0    no errors detected
!                          = 1    zero pivot encountered
!                                 (unsuccessful factorization)
!                          = 2    negative pivot encountered
!                                 (successful factorization)
!
!  
!
  integer   jt(1), iwksp(1)
  dimension d(1), t(ndim,1), wksp(1)
  logical   propa
!
!
  n = nn
  maxt = maxtt
  iflag = 0
  propa = ipropa == 1
  if (maxt < 1) go to 500
  nm1 = n - 1
  if (meth /= 1 .or. .not. propa) go to 20
!
!  ic, propa = t.
!
  do 15 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 10 j = 1,maxt
        kf = k + jt(j)
        if (kf <= n) d(kf) = d(kf) - t(k,j)**2/pivot
 10      continue
 15   continue
  if (d(n) == 0.0) go to 995
  go to 500
 20   if (meth /= 2 .or. .not. propa) go to 50
!
!  mic, propa = t.
!
  do 25 i = 1,n
 25   wksp(i) = 0.0
  do 35 j = 1,maxt
     do 30 i = 1,n
 30      wksp(i) = wksp(i) + t(i,j)
 35   continue
  do 45 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 40 i = 1,maxt
        kf = k + jt(i)
        if (kf > n) go to 40
        term = t(k,i)/pivot
        d(kf) = d(kf) - term*(omega*wksp(k)-(omega-1.0)*t(k,i))
 40      continue
 45   continue
  if (d(n) == 0.0) go to 995
  go to 500
!
!  ic, mic for propa = f.
!
 50   nbig = maxt + 1
  do 70 i = 1,maxt
     do 65 j = i,maxt
        if (j == i) go to 65
        iloc = (j - 1)*maxt + i
        id = iabs (jt(j) - jt(i))
        do 60 k = 1,maxt
           if (jt(k) /= id) go to 60
           iwksp(iloc) = k
           go to 65
 60         continue
        iwksp(iloc) = nbig
 65      continue
 70   continue
  do 100 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 95 i = 1,maxt
        kf = k + jt(i)
        if (kf > n) go to 95
        do 75 j = i,maxt
 75         wksp(j) = t(k,i)*t(k,j)/pivot
        d(kf) = d(kf) - wksp(i)
        do 90 j = i,maxt
           if (j == i) go to 90
           kg = k + jt(j)
           if (kg > n) go to 90
           iloc = (j-1)*maxt+i
           id = iwksp(iloc)
           if (id == nbig) go to 85
           kff = min (kf,kg)
           t(kff,id) = t(kff,id) - wksp(j)
           go to 90
 85            if (meth == 1) go to 90
           d(kf) = d(kf) - omega*wksp(j)
           d(kg) = d(kg) - omega*wksp(j)
 90         continue
 95      continue
 100  continue
  if (d(n) == 0.0) go to 995
!
!  store reciprocals of pivots.
!
 500  do 505 i = 1,n
 505  d(i) = 1.0/d(i)
  if (maxt < 1 .or. propa) go to 990
  do 515 j = 1,maxt
     len = n - jt(j)
     do 510 i = 1,len
 510     t(i,j) = d(i)*t(i,j)
 515  continue
!
!  check for negative pivots.
!
 990  if (vmin(n,d) < 0.0) iflag = 2
  return
!
!  error - matrix cannot be factored since a pivot is zero
!
 995  iflag = 1
  return
end
subroutine icfcp (ndimr,ndimi,nn,maxcc,jc,d,c,ncolor,nt,nb,meth,ipropa, &
  ipt,omega,iflag)
!
!*******************************************************************************
!
!! ICFCP computes an incomplete factorization.  (Purdue storage, multicolor)
!
!
!   The matrix is stored in d and c and the factorization replaces it.
!     
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          maxc   number of columns in c array
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          meth   point factorization wanted
!                  = 1   ic
!                  = 2   mic
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          ipt    integer pointer vector of length ncolor + 1
!          omega  modification factor between 0.0 and 1.0
!                  = 0     no modification
!                  = 1     full modification
!          iflag  indicator of factorization stability
!                    iflag = 0    no errors detected
!                          = 1    zero pivot encountered
!                                 (unsuccessful factorization)
!                          = 2    negative pivot encountered
!                                 (successful factorization)
!
!  
!
  integer   jc(ndimi,1), nt(1), nb(1), ipt(1)
  dimension d(1), c(ndimr,1)
  logical   propa
!
!
  n = nn
  maxc = maxcc
  ncol = ncolor
  iflag = 0
  propa = ipropa == 1
  if (maxc < 1) go to 75
!
!  do factorization.
!
  do 65 icol = 1,ncol-1
     k1 = ipt(icol) + 1
     k2 = ipt(icol+1)
     j22 = nt(icol)
     if (j22 <= 0) go to 65
     do 60 k = k1,k2
        pivot = d(k)
        if (pivot == 0.0) go to 995
        do 55 l1 = icol+1,ncol
           i1 = ipt(l1) + 1
           i2 = ipt(l1+1)
           j11 = nt(l1) + 1
           j12 = nt(l1) + nb(l1)
           j32 = nt(l1)
           if (j11 > j12) go to 55
           do 50 j1 = j11,j12
           do 45 i = i1,i2
              jcol1 = jc(i,j1)
              if (jcol1 /= k) go to 45
              term1 = c(i,j1)/pivot
              do 40 j2 = 1,j22
                 j = jc(k,j2)
                 if (j <= k) go to 40
                 term2 = term1*c(k,j2)
                 if (j == i) go to 35
                 if (propa) go to 30
                 if (j > i) go to 20
                 do 15 j3 = j11,j12
                    if (jc(i,j3) /= j) go to 15
                    c(i,j3) = c(i,j3) - term2
                    go to 40
 15                  continue
                 go to 30
 20                  if (j32 <= 0) go to 30
                 do 25 j3 = 1,j32
                    if (jc(i,j3) /= j) go to 25
                    c(i,j3) = c(i,j3) - term2
                    go to 40
 25                  continue
 30                  if (meth == 1) go to 40
 35                  d(i) = d(i) - omega*term2
 40               continue
 45            continue
 50            continue
 55         continue
 60      continue
 65   continue
  k1 = ipt(ncol) + 1
  k2 = ipt(ncol+1)
  do 70 k = k1,k2
     if (d(k) == 0.0) go to 995
 70   continue
!
!  store reciprocals of pivots.
!
 75   do 80 i = 1,n
 80   d(i) = 1.0/d(i)
  if (maxc < 1 .or. propa) go to 990
  do 105 icol = 1,ncol
     nt2 = nt(icol)
     i1 = ipt(icol) + 1
     i2 = ipt(icol+1)
     do 100 j = 1,maxc
        if (j > nt2) go to 90
        do 85 i = i1,i2
 85         c(i,j) = d(i)*c(i,j)
        go to 100
 90         do 95 i = i1,i2
 95         c(i,j) = c(i,j)*d(jc(i,j))
 100     continue
 105  continue
!
!  check for negative pivots.
!
 990  if (vmin(n,d) < 0.0) iflag = 2
  return
!
!  error - matrix cannot be factored since a pivot is zero
!
 995  iflag = 1
  return
end
subroutine icfn (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,meth, ipropa,omega,wksp, &
  iwksp,iflag)
!
!*******************************************************************************
!
!! ICFN computes an incomplete factorization.  (nonsymmetric diagonal storage)
!
!
!  The matrix is stored in d, t, and b and the factorization replaces it.
!     
!
!  Parameters:
!
!          ndim   row dimension of t,b arrays
!          n      order of system (= nn)
!          maxt   number of columns in t array
!          maxb   number of columns in b array
!          jt     integer vector giving the diagonal indices of
!                  the corresponding columns in t
!          jb     integer vector giving the diagonal indices of
!                  the corresponding columns in b
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          t      array of active size n by maxt giving the
!                  super-diagonals of the matrix
!          b      array of active size n by maxb giving the
!                  sub-diagonals of the matrix
!          meth   point factorization wanted
!                  = 1   ic
!                  = 2   mic
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          omega  modification factor between 0.0 and 1.0
!                  = 0     no modification
!                  = 1     full modification
!          wksp   workspace vector of length n
!          iwksp  integer workspace of length maxb*maxt
!          iflag  indicator of factorization stability
!                    iflag = 0    no errors detected
!                          = 1    zero pivot encountered
!                                 (unsuccessful factorization)
!                          = 2    negative pivot encountered
!                                 (successful factorization)
!
!  
!
  integer   jt(1), jb(1), iwksp(1)
  dimension d(1), t(ndim,1), b(ndim,1), wksp(1)
  logical   propa
!
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  iflag = 0
  propa = ipropa == 1
  if (maxt < 1  .or.  maxb < 1) go to 500
  nm1 = n - 1
  if (meth /= 1 .or. .not. propa) go to 30
!
!  ic, propa = t.
!
  nval = 0
  do 15 j = 1,maxb
     i1 = -jb(j)
     do 10 i = 1,maxt
        i2 = jt(i)
        if (i1 /= i2) go to 10
        nval = nval + 1
        iwksp(3*nval-2) = j
        iwksp(3*nval-1) = i
        iwksp(3*nval) = i2
        go to 15
 10      continue
 15   continue
  if (nval == 0) go to 500
  do 25 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 20 j = 1,nval
        kf = k + iwksp(3*j)
        if (kf > n) go to 20
        i1 = iwksp(3*j-2)
        i2 = iwksp(3*j-1)
        d(kf) = d(kf) - b(kf,i1)*t(k,i2)/pivot
 20      continue
 25   continue
  if (d(n) == 0.0) go to 995
  go to 500
 30   if (meth /= 2 .or. .not. propa) go to 70
!
!  mic, propa = t.
!
  do 35 i = 1,n
 35   wksp(i) = 0.0
  do 45 j = 1,maxt
     do 40 i = 1,n
 40      wksp(i) = wksp(i) + t(i,j)
 45   continue
  do 55 i = 1,maxb
     i1 = -jb(i)
     do 50 j = 1,maxt
        i2 = jt(j)
        if (i1 /= i2) go to 50
        iwksp(i) = j
        go to 55
 50      continue
     iwksp(i) = 0
 55   continue
  do 65 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 60 i = 1,maxb
        kf = k - jb(i)
        if (kf > n) go to 60
        term = b(kf,i)/pivot
        t1 = 0.0
        i1 = iwksp(i)
        if (i1 /= 0) t1 = t(k,i1)
        d(kf) = d(kf) - term*(omega*wksp(k)-(omega-1.0)*t1)
 60      continue
 65   continue
  if (d(n) == 0.0) go to 995
  go to 500
!
!  ic, mic for propa = f.
!
 70   nbig = maxt + maxb
  do 105 i = 1,maxb
     do 100 j = 1,maxt
        iloc = (j - 1)*maxb + i
        id = jt(j) + jb(i)
        if (id) 75,85,90
 75         do 80 k = 1,maxb
           if (jb(k) /= id) go to 80
           iwksp(iloc) = -k
           go to 100
 80         continue
        iwksp(iloc) = nbig
        go to 100
 85         iwksp(iloc) = 0
        go to 100
 90         do 95 k = 1,maxt
           if (jt(k) /= id) go to 95
           iwksp(iloc) = k
           go to 100
 95         continue
        iwksp(iloc) = nbig
 100     continue
 105  continue
  do 140 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 135 i = 1,maxb
        kf = k - jb(i)
        if (kf > n) go to 135
        do 110 j = 1,maxt
 110        wksp(j) = b(kf,i)*t(k,j)/pivot
        do 130 j = 1,maxt
           iloc = (j-1)*maxb+i
           id = iwksp(iloc)
           if (id) 115,120,125
 115           mid = -id
           b(kf,mid) = b(kf,mid) - wksp(j)
           go to 130
 120           d(kf) = d(kf) - wksp(j)
           go to 130
 125           if (id /= nbig) t(kf,id) = t(kf,id) - wksp(j)
           if (id == nbig .and. meth == 2) then
             d(kf) = d(kf) - omega*wksp(j)
           end if
 130        continue
 135     continue
 140  continue
  if (d(n) == 0.0) go to 995
!
!  store reciprocals of pivots.
!
 500  do 505 i = 1,n
 505  d(i) = 1.0/d(i)
  if (maxt < 1 .or. propa) go to 520
  do 515 j = 1,maxt
     len = n - jt(j)
     do 510 i = 1,len
 510     t(i,j) = d(i)*t(i,j)
 515  continue
 520  if (maxb < 1 .or. propa) go to 990
  do 530 j = 1,maxb
     ind = jb(j)
     len = n + ind
     do 525 i = 1,len
 525     b(i-ind,j) = d(i)*b(i-ind,j)
 530  continue
!
!  check for negative pivots.
!
 990  if (vmin(n,d) < 0.0) iflag = 2
  return
!
!  error - matrix cannot be factored since a pivot is zero
!
 995  iflag = 1
  return
end
subroutine icfnp (ndimr,ndimi,nn,maxtt,maxbb,jt,jb,d,t,b,meth,ipropa, &
  omega,iflag)
!
!*******************************************************************************
!
!! ICFNP computes an incomplete factorization.  (Purdue storage, nonsymmetric matrix)
!
!
!  The matrix is stored in d, t, and b and the factorization replaces it.
!     
!
!  Parameters:
!
!          ndimr  row dimension of t and b arrays
!          ndimi  row dimension of jt and jb arrays
!          n      order of system (= nn)
!          maxt   number of columns in t,jt arrays
!          maxb   number of columns in b,jb arrays
!          jt     integer array giving the column indices of the
!                  corresponding elements in t
!          jb     integer array giving the column indices of the
!                  corresponding elements in b
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          t      array of active size n by maxt giving the
!                  upper triangle of the matrix
!          b      array of active size n by maxb giving the
!                  lower triangle of the matrix
!          meth   point factorization wanted
!                  = 1   ic
!                  = 2   mic
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          omega  modification factor between 0.0 and 1.0
!                  = 0     no modification
!                  = 1     full modification
!          iflag  indicator of factorization stability
!                    iflag = 0    no errors detected
!                          = 1    zero pivot encountered
!                                 (unsuccessful factorization)
!                          = 2    negative pivot encountered
!                                 (successful factorization)
!
!  
!
  integer   jt(ndimi,1), jb(ndimi,1)
  dimension d(1), t(ndimr,1), b(ndimr,1)
  logical   propa
!
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  iflag = 0
  propa = ipropa == 1
!
  if (maxt < 1  .or.  maxb < 1) go to 50
  nm1 = n - 1
  do 45 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     kp1 = k + 1
     do 40 j1 = 1,maxb
     do 35 i = kp1,n
        jcol1 = jb(i,j1)
        if (jcol1 /= k) go to 35
        term1 = b(i,j1)/pivot
        do 30 j2 = 1,maxt
           j = jt(k,j2)
           if (j <= k) go to 30
           term2 = term1*t(k,j2)
           jdiff = j - i
           if (jdiff == 0) go to 27
           if (propa) go to 25
           if (jdiff > 0) go to 15
           do 10 j3 = 1,maxb
              if (jb(i,j3) /= j) go to 10
              b(i,j3) = b(i,j3) - term2
              go to 30
 10            continue
           go to 25
 15            do 20 j3 = 1,maxt
              if (jt(i,j3) /= j) go to 20
              t(i,j3) = t(i,j3) - term2
              go to 30
 20            continue
 25            if (meth == 1) go to 30
 27            d(i) = d(i) - omega*term2
 30         continue
 35      continue
 40      continue
 45   continue
  if (d(n) == 0.0) go to 995
!
!  store reciprocals of pivots.
!
 50   do 55 i = 1,n
 55   d(i) = 1.0/d(i)
  if (maxt < 1 .or. propa) go to 70
  do 65 j = 1,maxt
     do 60 i = 1,n
 60      t(i,j) = d(i)*t(i,j)
 65   continue
 70   if (maxb < 1 .or. propa) go to 990
  do 80 j = 1,maxb
     do 75 i = 1,n
 75      b(i,j) = b(i,j)*d(jb(i,j))
 80   continue
!
!  check for negative pivots.
!
 990  if (vmin(n,d) < 0.0) iflag = 2
  return
!
!  error - matrix cannot be factored since a pivot is zero
!
 995  iflag = 1
  return
end
subroutine icfp (ndimr,ndimi,nn,maxtt,jt,d,t,meth,ipropa,omega,wksp,iflag)
!
!*******************************************************************************
!
!! ICFP computes an incomplete factorization.  (Purdue storage, symmetric matrix)
!
!
!  The matrix is stored in d and t and the factorization replaces it.
!     
!
!  Parameters:
!
!          ndimr  row dimension of t array
!          ndimi  row dimension of jt array
!          n      order of system (= nn)
!          maxt   number of columns in t array
!          jt     integer array of active size n by maxt giving the
!                  column numbers of the corresponding elements in t
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          t      array of active size n by maxt giving the
!                  coefficients of the upper triangle of the matrix
!          meth   point factorization wanted
!                  = 1   ic
!                  = 2   mic
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          omega  modification factor between 0.0 and 1.0
!                  = 0     no modification
!                  = 1     full modification
!          wksp   workspace array of length n
!          iflag  indicator of factorization stability
!                    iflag = 0    no errors detected
!                          = 1    zero pivot encountered
!                                 (unsuccessful factorization)
!                          = 2    negative pivot encountered
!                                 (successful factorization)
!
!  
!
  dimension d(1), t(ndimr,1), wksp(1)
  integer   jt(ndimi,1)
  logical   propa
!
!
  n = nn
  maxt = maxtt
  iflag = 0
  propa = ipropa == 1
  if (maxt < 1) go to 500
  nm1 = n - 1
  if (meth /= 1 .or. .not. propa) go to 20
!
!  ic, propa = t.
!
  do 15 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 10 j = 1,maxt
        jcol = jt(k,j)
        d(jcol) = d(jcol) - t(k,j)**2/pivot
 10      continue
 15   continue
  if (d(n) == 0.0) go to 995
  go to 500
 20   if (meth /= 2 .or. .not. propa) go to 50
!
!  mic, propa = t.
!
  do 25 i = 1,n
 25   wksp(i) = 0.0
  do 35 j = 1,maxt
     do 30 i = 1,n
 30      wksp(i) = wksp(i) + t(i,j)
 35   continue
  do 45 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 40 i = 1,maxt
        jcol = jt(k,i)
        if (jcol == k) go to 40
        term = t(k,i)/pivot
        d(jcol) = d(jcol) - term*(omega*wksp(k)-(omega-1.0)*t(k,i))
 40      continue
 45   continue
  if (d(n) == 0.0) go to 995
  go to 500
!
!  ic, mic for propa = f.
!
 50   do 70 k = 1,nm1
     pivot = d(k)
     if (pivot == 0.0) go to 995
     do 65 j1 = 1,maxt
        jcol1 = jt(k,j1)
        if (jcol1 == k) go to 65
        d(jcol1) = d(jcol1) - (t(k,j1)**2)/pivot
        term1 = t(k,j1)/pivot
        do 60 j2 = 1,maxt
           jcol2 = jt(k,j2)
           if (jcol2 <= jcol1) go to 60
           if (jcol2 == k) go to 60
           term2 = term1*t(k,j2)
           do 55 j3 = 1,maxt
              if (jcol2 /= jt(jcol1,j3)) go to 55
              t(jcol1,j3) = t(jcol1,j3) - term2
              go to 60
 55            continue
           if (meth == 1) go to 60
           d(jcol1) = d(jcol1) - omega*term2
           d(jcol2) = d(jcol2) - omega*term2
 60         continue
 65      continue
 70   continue
  if (d(n) == 0.0) go to 995
!
!  store reciprocals of pivots and scale t.
!
 500  do 510 i = 1,n
 510  d(i) = 1.0/d(i)
  if (maxt < 1 .or. propa) go to 990
  do 520 j = 1,maxt
     do 515 i = 1,n
 515     t(i,j) = d(i)*t(i,j)
 520   continue
!
!  check for negative pivots.
!
 990  if (vmin(n,d) < 0.0) iflag = 2
  return
!
!  error - matrix cannot be factored since a pivot is zero
!
 995  iflag = 1
  return
end
subroutine icfs (ndim,nn,maxbb,jb,d,b,ipropa,irwise,iwksp,x)
!
!*******************************************************************************
!
!! ICFS does an IC forward solve (natural ordering, diagonal storage).
!
!
!        (i + b)*x = y    if not property a
!        (i + b*d)*x = y  if property a
!
!  Parameters:
!
!        ndim   row dimension of b array
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxb
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
  logical   propa
!
  n = nn
  maxb = maxbb
  propa = ipropa == 1
  if (maxb < 1) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 70
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxb
 15   iwksp(i) = 1 - jb(i)
!
!  determine nc, imin.
!
 20   nc = n
  do 25 i = 1,maxb
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 25
     nc = nterm
     imin = i
 25   continue
  if (nc >= n) return
  ndel = -jb(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 50
!
!  special case for first minor subdiagonal.
!
  nc1 = n
  do 30 i = 1,maxb
     if (i == imin) go to 30
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imin) = nc1 + 1
  if (propa) go to 40
  do 35 j = ibeg,nc1
 35   x(j) = x(j) - b(j,imin)*x(j-1)
  go to 20
 40   do 45 j = ibeg,nc1
 45   x(j) = x(j) - d(j-1)*b(j,imin)*x(j-1)
  go to 20
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imin) = iwksp(imin) + ndel
  iend = min (ibeg+ndel-1,n)
  if (propa) go to 60
!dir$ ivdep
  do 55 i = ibeg,iend
 55   x(i) = x(i) - b(i,imin)*x(i-ndel)
  go to 20
!dir$ ivdep
 60   do 65 i = ibeg,iend
 65   x(i) = x(i) - d(i-ndel)*b(i,imin)*x(i-ndel)
  go to 20
!
!  rowwise algorithm.
!
 70   if (propa) go to 90
  do 85 i = 2,n
     do 75 j = 1,maxb
 75      iwksp(j) = max (1,i+jb(j))
     sum = x(i)
     do 80 j = 1,maxb
 80      sum = sum - b(i,j)*x(iwksp(j))
     x(i) = sum
 85   continue
  return
 90   do 105 i = 2,n
     do 95 j = 1,maxb
 95      iwksp(j) = max (1,i+jb(j))
     sum = x(i)
     do 100 j = 1,maxb
 100     sum = sum - d(iwksp(j))*b(i,j)*x(iwksp(j))
     x(i) = sum
 105  continue
  return
end
subroutine icfscp (ndimr,ndimi,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,x)
!
!*******************************************************************************
!
!! ICFSCP does a forward IC solve.  (Purdue storage, multicolor)
!
!
!       (i + b)*x = y    if ipropa = 0
!       (d + b)*x = y    if ipropa = 1
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!          x      on input, x contains y
!                 on output, x is the solution to the back solve
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimr,1), x(1), wksp(1)
  logical   propa
!
  propa = ipropa == 1
!
  ist = 1
  do 25 icol = 1,ncolor
     npt = nc(icol)
     ied = ist + npt - 1
     j1 = nt(icol) + 1
     mj = nb(icol)
     call vsubp (ndimr,ndimi,npt,mj,c(ist,j1),jc(ist,j1),x(ist),x,wksp)
     if (.not. propa) go to 20
     do 15 i = ist,ied
 15      x(i) = x(i)*d(i)
 20      ist = ist + npt
 25   continue
  return
end
subroutine icfsct (ndimr,ndimi,jc,d,c,ncolor,nc,nt,ipropa,wksp,x)
!
!*******************************************************************************
!
!! ICFSCT does a transpose forward ic solve.  (Purdue storage, multicolor)
!
!
!     (i + (t**t))*x = y    if ipropa = 0
!     (d + (t**t))*x = y    if ipropa = 1
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length max(nc(i))
!          x      on input, x contains y
!                 on output, x is the solution to the forward solve
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1)
  dimension d(1), c(ndimr,1), x(1), wksp(1)
  logical   propa
!
  propa = ipropa == 1
!
  ist =  1
  do 25 icol = 1,ncolor
     npt = nc(icol)
     ied = ist + npt - 1
     if (.not. propa) go to 20
     do 15 i = ist,ied
 15      x(i) = x(i)*d(i)
 20      j2 = nt(icol)
     call vsubpt (ndimr,ndimi,npt,j2,c(ist,1),jc(ist,1),x,x(ist),wksp)
     ist = ist + npt
 25   continue
  return
end
subroutine icfsp (ndimr,ndimi,n,maxb,jb,d,b,ipropa,x)
!
!*******************************************************************************
!
!! ICFSP does an IC forward solve (natural ordering, Purdue storage).
!
!
!        (i + b)*x = y    if ipropa = 0
!        (d + b)*x = y    if ipropa = 1
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the lower
!                triangle of the factorization if ipropa = 0
!                or the lower triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), b(ndimr,1)
  integer   jb(ndimi,1)
  logical   propa
!
  propa = ipropa == 1
  if (maxb >= 1) go to 15
  if (.not. propa) return
  do 10 i = 1,n
 10   x(i) = x(i)*d(i)
  return
 15   do 25 i = 1,n
     sum = x(i)
     do 20 j = 1,maxb
        sum = sum - b(i,j)*x(jb(i,j))
 20      continue
     if (propa) sum = sum*d(i)
     x(i) = sum
 25   continue
  return
end
subroutine icfst (ndim,nn,maxtt,jt,d,t,ipropa,irwise,iwksp,x)
!
!*******************************************************************************
!
!! ICFST does an IC forward solve (natural ordering, diagonal storage).
!
!
!        (i + (t**t))*x = y    if not property a
!        (i + (t**t)*d)*x = y  if property a
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
  logical   propa
!
  n = nn
  maxt = maxtt
  nm1 = n - 1
  propa = ipropa == 1
  if (maxt < 1) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 70
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxt
 15   iwksp(i) = jt(i) + 1
!
!  determine nc, imin.
!
 20   nc = n
  do 25 i = 1,maxt
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 25
     nc = nterm
     imin = i
 25   continue
  if (nc >= n) return
  ndel = jt(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 50
!
!  special case for first minor subdiagonal.
!
  nc1 = n
  do 30 i = 1,maxt
     if (i == imin) go to 30
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imin) = nc1 + 1
  if (propa) go to 40
  do 35 j = ibeg,nc1
 35   x(j) = x(j) - t(j-1,imin)*x(j-1)
  go to 20
 40   do 45 j = ibeg,nc1
 45   x(j) = x(j) - d(j-1)*t(j-1,imin)*x(j-1)
  go to 20
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imin) = iwksp(imin) + ndel
  iend = min (ibeg+ndel-1,n)
  if (propa) go to 60
!dir$ ivdep
  do 55 i = ibeg,iend
 55   x(i) = x(i) - t(i-ndel,imin)*x(i-ndel)
  go to 20
!dir$ ivdep
 60   do 65 i = ibeg,iend
 65   x(i) = x(i) - d(i-ndel)*t(i-ndel,imin)*x(i-ndel)
  go to 20
!
!  rowwise algorithm.
!
 70   do 85 i = 1,nm1
     do 75 j = 1,maxt
 75      iwksp(j) = min (n,i+jt(j))
     term = x(i)
     if (propa) term = term*d(i)
     do 80 j = 1,maxt
 80      x(iwksp(j)) = x(iwksp(j)) - t(i,j)*term
 85   continue
  return
end
subroutine icfstp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
!
!*******************************************************************************
!
!! ICFSTP does a transpose IC forward solve (natural ordering, Purdue storage).
!
!
!        (i + (t**t))*x = y    if ipropa = 0
!        (d + (t**t))*x = y    if ipropa = 1
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
  logical   propa
!
  propa = ipropa == 1
  if (maxt >= 1) go to 15
  if (.not. propa) return
  do 10 i = 1,n
 10   x(i) = x(i)*d(i)
  return
 15   do 25 i = 1,n
     if (propa) x(i) = x(i)*d(i)
     term = x(i)
     do 20 j = 1,maxt
        x(jt(i,j)) = x(jt(i,j)) - t(i,j)*term
 20      continue
 25   continue
  return
end
subroutine icfv (ndim,nn,maxtt,jt,d,t,meth,ipropa,omega,wksp,iwksp,iflag)
!
!*******************************************************************************
!
!! ICFV computes an incomplete factorization.  (symmetric diagonal storage, vectorized version)
!
!
!  The matrix is stored in d and t and the factorization replaces it.
!     
!
!  Parameters:
!
!          ndim   row dimension of t array
!          n      order of system (= nn)
!          maxt   number of columns in t array
!          jt     integer vector giving the diagonal indices of
!                  the corresponding columns in t
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          t      array of active size n by maxt giving the
!                  super-diagonals of the matrix
!          meth   point factorization wanted
!                  = 1   ic
!                  = 2   mic
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          omega  modification factor between 0.0 and 1.0
!                  = 0     no modification
!                  = 1     full modification
!          wksp   workspace vector of length n
!          iwksp  integer workspace of length maxt**2
!          iflag  indicator of factorization stability
!                    iflag = 0    no errors detected
!                          = 1    zero pivot encountered
!                                 (unsuccessful factorization)
!                          = 2    negative pivot encountered
!                                 (successful factorization)
!
!  
!
  integer   jt(1), iwksp(1)
  dimension d(1), t(ndim,1), wksp(1)
  logical propa
!
!
  n = nn
  maxt = maxtt
  iflag = 0
  propa = ipropa == 1
  if (maxt < 1) go to 500
  if (meth /= 1 .or. .not. propa) go to 45
!
!  ic, propa = t.
!
  do 10 i = 1,maxt
 10   iwksp(i) = jt(i) + 1
!
!  determine nc, imin.
!
 15   nc = n
  do 20 i = 1,maxt
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 20
     nc = nterm
     imin = i
 20   continue
  if (nc >= n) go to 500
  ndel = jt(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 35
!
!  special case for first super-diagonal.
!
  nc1 = n
  do 25 i = 1,maxt
     if (i == imin) go to 25
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 25   continue
  iwksp(imin) = nc1 + 1
  do 30 j = ibeg,nc1
 30   d(j) = d(j) - (t(j-1,imin)**2)/d(j-1)
  go to 15
!
!  far diagonals.
!
 35   iwksp(imin) = iwksp(imin) + ndel
  ied = min (ibeg+ndel-1,n)
!dir$ ivdep
  do 40 i = ibeg,ied
 40   d(i) = d(i) - (t(i-ndel,imin)**2)/d(i-ndel)
  go to 15
 45   if (meth /= 2 .or. .not. propa) go to 100
!
!  mic, propa = t.
!
  do 50 i = 1,n
 50   wksp(i) = 0.0
  do 60 j = 1,maxt
     do 55 i = 1,n
 55      wksp(i) = wksp(i) + t(i,j)
 60   continue
  do 65 i = 1,maxt
 65   iwksp(i) = jt(i) + 1
!
!  determine nc, imin.
!
 70   nc = n
  do 75 i = 1,maxt
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 75
     nc = nterm
     imin = i
 75   continue
  if (nc >= n) go to 500
  ndel = jt(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 90
!
!  special case for first super-diagonal.
!
  nc1 = n
  do 80 i = 1,maxt
     if (i == imin) go to 80
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 80   continue
  iwksp(imin) = nc1 + 1
  do 85 j = ibeg,nc1
 85   d(j) = d(j) - t(j-1,imin)*(omega*wksp(j-1)-(omega-1.0)*t(j-1,imin))/d(j-1)
  go to 70
!
!  far diagonals.
!
 90   iwksp(imin) = iwksp(imin) + ndel
  ied = min (ibeg+ndel-1,n)
!dir$ ivdep
  do 95 i = ibeg,ied
 95   d(i) = d(i) - t(i-ndel,imin)*(omega*wksp(i-ndel)- &
        (omega-1.0)*t(i-ndel,imin))/d(i-ndel)
  go to 70
!
!  set up pointers for propa = f case.
!
 100  nbig = maxt + 1
  do 115 i = 1,maxt
     do 110 j = 1,maxt
        iloc = j*maxt + i
        id = iabs (jt(j) - jt(i))
        do 105 k = 1,maxt
           if (jt(k) /= id) go to 105
           iwksp(iloc) = k
           go to 110
 105        continue
        iwksp(iloc) = nbig
 110     continue
 115  continue
!
!  ic, mic for propa = f.
!
  do 120 i = 1,maxt
 120  iwksp(i) = jt(i) + 1
!
!  determine nc, imin.
!
 125  nc = n
  do 130 i = 1,maxt
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 130
     nc = nterm
     imin = i
 130  continue
  if (nc >= n) go to 500
  ndel = jt(imin)
  iwksp(imin) = iwksp(imin) + ndel
  ibeg = nc + 1
  ied = min (ibeg+ndel-1,n)
!dir$ ivdep
  do 135 i = ibeg,ied
 135  d(i) = d(i) - (t(i-ndel,imin)**2)/d(i-ndel)
  do 160 j = 1,maxt
     jcol = jt(j)
     if (jcol <= ndel) go to 160
     iloc = j*maxt + imin
     id = iwksp(iloc)
     ied1 = min (ied,n-jcol+ndel)
     if (id == nbig) go to 145
!dir$ ivdep
     do 140 i = ibeg,ied1
 140     t(i,id) = t(i,id) - t(i-ndel,imin)*t(i-ndel,j)/d(i-ndel)
     go to 160
 145     if (meth == 1) go to 160
     do 150 i = ibeg,ied1
 150     wksp(i) = omega*t(i-ndel,imin)*t(i-ndel,j)/d(i-ndel)
     ish = jcol - ndel
     do 155 i = ibeg,ied1
        d(i) = d(i) - wksp(i)
        d(i+ish) = d(i+ish) - wksp(i)
 155     continue
 160  continue
  go to 125
!
!  store reciprocals of pivots.
!
 500  do 505 i = 1,n
     if (d(i) == 0.0) go to 995
 505  continue
  do 510 i = 1,n
 510  d(i) = 1.0/d(i)
  if (maxt < 1 .or. propa) go to 990
  do 520 j = 1,maxt
     len = n - jt(j)
     do 515 i = 1,len
 515     t(i,j) = d(i)*t(i,j)
 520  continue
!
!  check for negative pivots.
!
 990  if (vmin(n,d) < 0.0) iflag = 2
  return
!
!  error - matrix cannot be factored since a pivot is zero
!
 995  iflag = 1
  return
end
subroutine ics (ndim,nn,maxtt,jt,d,t,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICS does an IC solution (natural ordering, symmetric diagonal storage).
!
!
!        (i + (t**t))*inv(d)*(i + t)*x = y            propa = .false.
!        (i + (t**t)*d)*inv(d)*(i + d*t)*x = y        propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfst (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call icbs (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  return
end
subroutine ics1 (ndim,nn,maxtt,jt,d,t,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICS1 does an IC forward solution (natural ordering, symmetric diagonal storage).
!
!
!        (i + (t**t))*inv(d)*x = y            propa = .false.
!        (i + (t**t)*d)*inv(d)*x = y          propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfst (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = sqrt(abs(d(i)))*x(i)
  return
end
subroutine ics2 (ndim,nn,maxtt,jt,d,t,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICS2 does an IC back solution (natural ordering, symmetric diagonal storage).
!
!
!        (i + t)*x = y            propa = .false.
!        (i + d*t)*x = y          propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  do 10 i = 1,n
 10   x(i) = y(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  call icbs (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  return
end
subroutine ics3 (ndim,nn,maxtt,jt,d,t,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICS3 does an IC transpose backward solution (natural ordering, symmetric diagonal storage).
!
!
!        inv(d)*(i + t)*x = y                 propa = .false.
!        inv(d)*(i + d*t)*x = y               propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  do 15 i = 1,n
 15   x(i) = sqrt(abs(d(i)))*y(i)
  call icbs (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  return
end
subroutine ics4 (ndim,nn,maxtt,jt,d,t,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICS4 does an IC transpose forward solution (natural ordering, symmetric diagonal storage).
!
!
!        (i + (t**t))*x = y            propa = .false.
!        (i + (t**t)*d)*x = y          propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfst (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = x(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  return
end
subroutine icscp (ndimr,ndimi,nn,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,y,x)
!
!*******************************************************************************
!
!! ICSCP does an IC solve. (Purdue storage, multicolor)
!
!
!      (i + b)*d*(i + t)*x = y         if ipropa = 0
!      (d + b)*inv(d)*(d + t)*x = y    if ipropa = 1
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!          y      on input, y is the right-hand-side vector
!          x      on output, x is the solution to the forward solve
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimr,1), x(1), y(1), wksp(1)
!
  n = nn
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfscp (ndimr,ndimi,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*d(i)
  go to 30
 20   do 25 i = 1,n
 25   x(i) = x(i)/d(i)
 30   continue
  call icbscp (ndimr,ndimi,n,jc,d,c,ncolor,nc,nt,ipropa,wksp,x)
  return
end
subroutine icscp1 (ndimr,ndimi,nn,jc,d,c,ncolor,nc,nt,nb,ipropa, wksp,y,x)
!
!*******************************************************************************
!
!! ICSCP1 does an IC forward solve.  (Purdue storage, multicolor)
!
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!          y      on input, y is the right-hand-side vector
!          x      on output, x is the solution to the forward solve
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimr,1), x(1), y(1), wksp(1)
!
  n = nn
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfscp (ndimr,ndimi,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*sqrt(abs(d(i)))
  return
 20   do 25 i = 1,n
 25   x(i) = x(i)/sqrt(abs(d(i)))
  return
end
subroutine icscp2 (ndimr,ndimi,nn,jc,d,c,ncolor,nc,nt,ipropa,wksp,y,x)
!
!*******************************************************************************
!
!! ICSCP2 does an IC back solve.  (Purdue storage, multicolor)
!
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!          y      on input, y is the right-hand-side vector
!          x      on output, x is the solution to the forward solve
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1)
  dimension d(1), c(ndimr,1), x(1), y(1), wksp(1)
!
  n = nn
!
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = y(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  go to 30
 20   do 25 i = 1,n
 25   x(i) = y(i)/(sign(1.0,d(i))*sqrt(abs(d(i))))
 30   continue
  call icbscp (ndimr,ndimi,n,jc,d,c,ncolor,nc,nt,ipropa,wksp,x)
  return
end
subroutine icscp3 (ndimr,ndimi,nn,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,y,x)
!
!*******************************************************************************
!
!! ICSCP3 does a transpose IC forward solve.  (Purdue storage, multicolor)
!
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length max(nc(i))
!          y      on input, y is the right-hand-side vector
!          x      on output, x is the solution vector
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimr,1), x(1), y(1), wksp(1)
!
  n = nn
!
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = y(i)*sqrt(abs(d(i)))
  go to 30
 20   do 25 i = 1,n
 25   x(i) = y(i)/sqrt(abs(d(i)))
 30   continue
  call icbsct (ndimr,ndimi,n,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,x)
  return
end
subroutine icscp4 (ndimr,ndimi,nn,jc,d,c,ncolor,nc,nt,ipropa,wksp,y,x)
!
!*******************************************************************************
!
!! ICSCP4 does a transpose IC back solve.  (Purdue storage, multicolor)
!
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length max(nc(i))
!          y      on input, y is the right-hand-side vector
!          x      on output, x is the solution vector
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1)
  dimension d(1), c(ndimr,1), x(1), y(1), wksp(1)
!
  n = nn
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfsct (ndimr,ndimi,jc,d,c,ncolor,nc,nt,ipropa,wksp,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  return
 20   do 25 i = 1,n
 25   x(i) = x(i)/(sign(1.0,d(i))*sqrt(abs(d(i))))
  return
end
subroutine icscpt (ndimr,ndimi,nn,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,y,x)
!
!*******************************************************************************
!
!! ICSCPT does a transpose IC solve.  (Purdue storage, multicolor)
!
!
!     (i + (t**t))*d*(i + (b**t))*x = y       if ipropa = 0
!     (d + (t**t))*inv(d)*(d + (b**t))*x = y  if ipropa = 1
!
!  Parameters:
!
!          ndimr  row dimension of c array
!          ndimi  row dimension of jc array
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          ipropa property a flag
!                  = 0   matrix does not have property a
!                  = 1   matrix has property a
!          wksp   workspace vector of length max(nc(i))
!          y      on input, y is the right-hand-side vector
!          x      on output, x is the solution vector
!
!  
!
  integer   jc(ndimi,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimr,1), x(1), y(1), wksp(1)
!
  n = nn
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfsct (ndimr,ndimi,jc,d,c,ncolor,nc,nt,ipropa,wksp,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*d(i)
  go to 30
 20   do 25 i = 1,n
 25   x(i) = x(i)/d(i)
 30   continue
  call icbsct (ndimr,ndimi,n,jc,d,c,ncolor,nc,nt,nb,ipropa,wksp,x)
  return
end
subroutine icsn (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICSN does an IC solution (natural ordering, nonsymmetric diagonal storage).
!
!
!        (i + b)*inv(d)*(i + t)*x = y            propa = .false.
!        (i + b*d)*inv(d)*(i + d*t)*x = y        propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        b      array of active size n by maxb giving the sub-
!                diagonals of the factorization if not property a
!                or the sub-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1), b(ndim,1)
  integer   jt(1), jb(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfs (ndim,n,maxb,jb,d,b,ipropa,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call icbs (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  return
end
subroutine icsn1 (ndim,n,maxb,jb,d,b,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICSN1 does an IC forward pass (natural ordering, nonsymmetric diagonal storage).
!
!
!        (i + b)*inv(d)*(i + t)*x = y            propa = .false.
!        (i + b*d)*inv(d)*(i + d*t)*x = y        propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the sub-
!                diagonals of the factorization if not property a
!                or the sub-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfs (ndim,n,maxb,jb,d,b,ipropa,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = sqrt(abs(d(i)))*x(i)
  return
end
subroutine icsn2 (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICSN2 does an IC back pass (natural ordering, nonsymmetric diagonal storage).
!
!
!        (i + b)*inv(d)*(i + t)*x = y            propa = .false.
!        (i + b*d)*inv(d)*(i + d*t)*x = y        propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  call icbs (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  return
end
subroutine icsn3 (ndim,n,maxb,jb,d,b,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICSN3 does an IC transpose back pass (natural ordering, nonsymmetric diagonal storage).
!
!
!        (i + b)*inv(d)*(i + t)*x = y            propa = .false.
!        (i + b*d)*inv(d)*(i + d*t)*x = y        propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the sub-
!                diagonals of the factorization if not property a
!                or the sub-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
!
  do 15 i = 1,n
 15   x(i) = sqrt(abs(d(i)))*y(i)
  call icbst (ndim,n,maxb,jb,d,b,ipropa,irwise,iwksp,x)
  return
end
subroutine icsn4 (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICSN4 does an IC transpose forward pass (natural ordering, nonsymmetric diagonal storage).
!
!
!        (i + b)*inv(d)*(i + t)*x = y            propa = .false.
!        (i + b*d)*inv(d)*(i + d*t)*x = y        propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfst (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = x(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  return
end
subroutine icsnp (ndimr,ndimi,nn,maxtt,maxbb,jt,jb,d,t,b,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSNP does an IC solution (natural ordering, Purdue storage, nonsymmetric matrix).
!
!
!        (i + b)*d*(i + t)*x = y                  if ipropa = 0
!        (d + b)*inv(d)*(d + t)*x = y             if ipropa = 1
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        b      array of active size n by maxb giving the lower
!                triangle of the factorization if ipropa = 0
!                or the lower triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1), b(ndimr,1)
  integer   jt(ndimi,1), jb(ndimi,1)
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfsp (ndimr,ndimi,n,maxb,jb,d,b,ipropa,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*d(i)
  go to 30
 20   do 25 i = 1,n
 25   x(i) = x(i)/d(i)
 30   continue
  call icbsp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  return
end
subroutine icsnp1 (ndimr,ndimi,nn,maxb,jb,d,b,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSNP1 does an IC forward solution (natural ordering, Purdue storage, nonsymmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the lower
!                triangle of the factorization if ipropa = 0
!                or the lower triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndimr,1)
  integer   jb(ndimi,1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfsp (ndimr,ndimi,n,maxb,jb,d,b,ipropa,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*sqrt(abs(d(i)))
  return
 20   do 25 i = 1,n
 25   x(i) = x(i)/sqrt(abs(d(i)))
  return
end
subroutine icsnp2 (ndimr,ndimi,n,maxt,jt,d,t,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSNP2 does an IC back solution (natural ordering, Purdue storage, nonsymmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
!
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = y(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  go to 30
 20   do 25 i = 1,n
 25   x(i) = y(i)/(sign(1.0,d(i))*sqrt(abs(d(i))))
 30   continue
  call icbsp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  return
end
subroutine icsnp3 (ndimr,ndimi,n,maxb,jb,d,b,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSNP3 does a transpose IC forward solution (natural ordering, Purdue storage, nonsymmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        b      array of active size n by maxb giving the lower
!                triangle of the factorization if ipropa = 0
!                or the lower triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndimr,1)
  integer   jb(ndimi,1)
!
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = y(i)*sqrt(abs(d(i)))
  go to 30
 20   do 25 i = 1,n
 25   x(i)  = y(i)/sqrt(abs(d(i)))
 30   continue
  call icbstp (ndimr,ndimi,n,maxb,jb,d,b,ipropa,x)
  return
end
subroutine icsnp4 (ndimr,ndimi,n,maxt,jt,d,t,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSNP4 does a transpose IC back solution (natural ordering, Purdue storage, nonsymmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfstp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  return
 20   do 25 i = 1,n
 25   x(i) = x(i)/(sign(1.0,d(i))*sqrt(abs(d(i))))
  return
end
subroutine icsnt (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,ipropa,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! ICSNT does a transpose IC solution (natural ordering, nonsymmetric diagonal storage).
!
!
!       (i + (t**t))*inv(d)*(i + (b**t))*x = y        propa = .false.
!       (i + (t**t)*d)*inv(d)*(i + d*(b**t))*x = y    propa = .true.
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the super-
!                diagonals of the factorization if not property a
!                or the super-diagonals of the matrix if property a
!        b      array of active size n by maxb giving the sub-
!                diagonals of the factorization if not property a
!                or the sub-diagonals of the matrix if property a
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1), b(ndim,1)
  integer   jt(1), jb(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfst (ndim,n,maxt,jt,d,t,ipropa,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call icbst (ndim,n,maxb,jb,d,b,ipropa,irwise,iwksp,x)
  return
end
subroutine icsntp (ndimr,ndimi,nn,maxtt,maxbb,jt,jb,d,t,b,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSNTP does a transpose IC solution (natural ordering, Purdue storage, nonsymmetric matrix).
!
!
!        (i + (t**t))*d*(i + (b**t))*x = y        if ipropa = 0
!        (d + (t**t))*inv(d)*(d + (b**t))*x = y   if ipropa = 1
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        b      array of active size n by maxb giving the lower
!                triangle of the factorization if ipropa = 0
!                or the lower triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1), b(ndimr,1)
  integer   jt(ndimi,1), jb(ndimi,1)
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfstp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*d(i)
  go to 30
 20   do 25 i = 1,n
 25   x(i) = x(i)/d(i)
 30   continue
  call icbstp (ndimr,ndimi,n,maxb,jb,d,b,ipropa,x)
  return
end
subroutine icsp (ndimr,ndimi,nn,maxtt,jt,d,t,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSP does an IC solution (natural ordering, Purdue storage, symmetric matrix).
!
!
!        (i + (t**t))*d*(i + t)*x = y             if ipropa = 0
!        (d + (t**t))*inv(d)*(d + t)*x = y        if ipropa = 1
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
!
  n = nn
  maxt = maxtt
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfstp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*d(i)
  go to 30
 20   do 25 i = 1,n
 25   x(i) = x(i)/d(i)
 30   continue
  call icbsp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  return
end
subroutine icsp1 (ndimr,ndimi,nn,maxt,jt,d,t,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSP1 does an IC forward solution (natural ordering, Purdue storage, symmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call icfstp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = x(i)*sqrt(abs(d(i)))
  return
 20   do 25 i = 1,n
 25   x(i) = x(i)/sqrt(abs(d(i)))
  return
end
subroutine icsp2 (ndimr,ndimi,n,maxt,jt,d,t,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSP2 does an IC back solution (natural ordering, Purdue storage, symmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
!
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = y(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  go to 30
 20   do 25 i = 1,n
 25   x(i) = y(i)/(sign(1.0,d(i))*sqrt(abs(d(i))))
 30   continue
  call icbsp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  return
end
subroutine icsp3 (ndimr,ndimi,n,maxt,jt,d,t,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSP3 does an IC transpose forward solution (natural ordering, Purdue storage, symmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
!
  if (ipropa == 1) go to 20
  do 15 i = 1,n
 15   x(i) = y(i)*sqrt(abs(d(i)))
  go to 30
 20   do 25 i = 1,n
 25   x(i) = y(i)/sqrt(abs(d(i)))
 30   continue
  call icbsp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  return
end
subroutine icsp4 (ndimr,ndimi,n,maxt,jt,d,t,ipropa,y,x)
!
!*******************************************************************************
!
!! ICSP4 does an IC transpose back solution (natural ordering, Purdue storage, symmetric matrix).
!
!
!  Parameters:
!
!        ndimr  row dimension of real arrays
!        ndimi  row dimension of integer arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the factorization
!        t      array of active size n by maxt giving the upper
!                triangle of the factorization if ipropa = 0
!                or the upper triangle of the matrix if ipropa = 1
!        ipropa property a switch
!                = 0  matrix does not have property a
!                = 1  matrix does have property a
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndimr,1)
  integer   jt(ndimi,1)
!
  do i = 1,n
    x(i) = y(i)
  end do

  call icfstp (ndimr,ndimi,n,maxt,jt,d,t,ipropa,x)
  if (ipropa == 1) go to 20

  do 15 i = 1,n
 15   x(i) = x(i)*sign(1.0,d(i))*sqrt(abs(d(i)))
  return

 20   do 25 i = 1,n
 25   x(i) = x(i)/(sign(1.0,d(i))*sqrt(abs(d(i))))

  return
end
subroutine inithv (icall)
!
!*******************************************************************************
!
!! INITHV initializes dot and vector "haves" to FALSE.
!
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
  udhav  = .false.
  rdhav  = .false.
  rzhav  = .false.
  rzthav = .false.
  zdhav  = .false.
  zzthav = .false.
  ztdhav = .false.
  if (icall == 1) return
  rhave  = .false.
  zhave  = .false.
  zthave = .false.
!
  return
end
subroutine iom (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! IOM is the user interface to the (truncated) IOM algorithm.  
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call iomw (suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs,wksp(irpnt), &
    nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine iomw (suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
  wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! IOMW runs the (truncated) IOM algorithm.  
!
!
!  the reference is
! youcef saad, "krylov subspace methods.", mathematics of
! computation, vol. 37, july 1981, pp. 105f.
!
! in the symmetric case this algorithm reduces to the symmlq
! algorithm of paige and saunders, except paige and saunders have
! implemented a trick to avoid breakdown before convergence.  this
! trick is not implemented here.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  integer idotw, vect1, vect2, dots1, dots2, os
  logical uneed
  external suba, subql, subqr
  dimension iparm(30), rparm(30)
  dimension gdum(1), wkxxx(1)
  logical iql, iqr
  logical exact, gamize
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! next, the indexing functions.
!
  indv1(i) = vect1 + mod(i,nv)*n
  indbe2(i) = ibeta2 + mod(i,os)
  indc(i) = icos + mod(i,os)
  inds(i) = isin + mod(i,os)
  indu(i) = iu + mod(i,os+1)
  indw(i) = iw + n*mod(i,os)
!
! preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 10
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 996
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  gamize = .true.
  if (iqr) go to 995
  if (level >= 2) write (nout,496)
496   format (' iom')
! the following flag tells us whether the truncating actually
! throws out important information.  it should actually be set to
! true if the matrix is symmetric.
  exact = .false.
!
! initialize the stopping test.
!
  call inithv (0)
  zdhav = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
!  associated integer variables.
!
  os = iabs(ns1)
  iv = 1
  nv = os
  idotw = 1
  iw = 1
  vect1 = iw + iv*n*os
  vect2 = vect1
  dots1 = vect2 + iv*n*nv
  dots2 = dots1
  ibeta1 = dots2 + idotw*os
  ibeta2 = ibeta1
  icos = ibeta2 + os
  isin = icos + os
  iu = isin + os
  iv1 = iu + os+1
  iv2 = iv1 + n
  nwusd = max(nwusd,iv2-1+n)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  uneed = rcalp .or. zcalp .or. ztcalp .or. udhav .or. ntest == 6 .or. &
    level >= 3
!
!  Begin iteration loop.
!
! perform first-iterate calculations.
!
 10   if (is /= 0) go to 100
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  call pvec (n,nv,iv,1,os,idotw,is,1,1,wk(vect1),wk(dots1),0,wk(ibeta1), &
    gdum,gamize,wk(iv2),wkxxx,ier)
  gamma1 = gdum(1)
  if (ier < 0) go to 997
  gamma2 = gamma1
  vnorm1 = 1.0/gamma1
  vnorm2 = 1.0/gamma2
  zdot = vnorm1**2
  ucnp1= 0.0
!
 100  call inithv (1)
  zdhav = .true.
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
!
!  compute q(n+1), etc -- the direction vectors
!
  call suba (coef,jcoef,wfac,jwfac,n,wk(indv1(is)),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  call pvec (n,nv,iv,1,os,idotw,is+1,1,1,wk(vect1),wk(dots1),0,wk(ibeta1), &
    gdum,gamize,wk(iv2),wkxxx,ier)
  gamma1 = gdum(1)
  if (ier < 0) go to 997
  gamma2 = gamma1
!
!  now record norms.
!
  vn1old = vnorm1
  vnorm1 = 1.0/gamma1
  vn2old = vnorm2
  vnorm2 = 1.0/gamma2
!
!  now update the factorization
!
  ucnbar = ucnp1
  ibgn = max(0,is+1-os)
  do 1 i = ibgn,is
 1    wk(indu(i+1)) = -wk(indbe2(i))
  if (ibgn > 0) wk(indu(ibgn))= 0.0
  call qrupd (is+1,os+1,os,wk(icos),wk(isin),ucnbar,ucn,wk(iu),vn2old,ier)
  if (ier < 0) go to 998
  ucnp1 = wk(indu(is+1))
!
!  update the old w vector.
!
  if (is /= 0) then
    call vtriad (n,wk(indw(is-1)),xxx,ucnbar/ucn,wk(indw(is-1)),2)
  end if
!
!  now generate the new w vector.
!
  if (abs(ucnp1) < srelpr) go to 998
  call vcopy (n,wk(indv1(is)),wk(iv1))
  ibgn = max(1,is-os+1)
  iend = is
  if (iend < ibgn) go to 200
  do 201 i = ibgn,iend
 201  call vtriad (n,wk(iv1),wk(iv1),-wk(indu(i)),wk(indw(i-1)),1)
 200  continue
  call vtriad (n,wk(indw(is)),xxx,1.0/ucnp1,wk(iv1),2)
  if (is /= 0) go to 205
!
!  update iterate u(0).
!
  zold= 0.0
  zbar = vn1old
  if (uneed) call vtriad (n,u,u,zbar,wk(indw(0)),1)
  go to 210
!
!  update subsequent iterates u(n).
!
 205  zold = wk(indc(is))*zbar
  zbold = zbar
  zbar =-wk(inds(is))*zbar
  factor = zold
  if (uneed) factor = factor - zbold*ucn/ucnbar
  call vtriad (n,u,u,factor,wk(indw(is-1)),1)
  if (uneed) call vtriad (n,u,u,zbar,wk(indw(is)),1)
! to avoid breakdown for the symmetric indefinite case, we'd only add
! in w(is-1) here, i believe.
 210  continue
  zdot = (zbar/ucnp1*vnorm1)**2
!
! proceed to next iteration
!
  in = in + 1
  is = is + 1
  go to 10
!
!  Finish up.
!
 900  if (.not. uneed) call vtriad (n,u,u,zbar,wk(indw(is-1)),1)
  if (halt) go to 715
  ier = 1
  call ershow (ier,'iomw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' iom converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 995  ier = -16
  call ershow (ier,'iomw')
  return

 996  call ershow (ier,'iomw')
  go to 735

 997  ier = -13
  call ershow (ier,'iomw')
  go to 725

 998  ier = -14
  call ershow (ier,'iomw')
  go to 725

 999  ier = -2
  call ershow (ier,'iomw')
  go to 735

end
function ipstr (omega)
!
!*******************************************************************************
!
!! IPSTR finds a suitable exponent for OMEGA-1.
!
!
!  Discussion:
!
!    IPSTR is the smallest integer such that
!
!      ipstr * (omega-1)**(ipstr-1) <= 0.50. 
!
!    IPSTR is required to be greater than 5.
!
!  Parameters:
!
!    omega  relaxation factor for sor method
!
  wm1 = omega - 1.0
  factor = wm1**5
!
  do 10 ip = 6,940
     if ( float (ip)*factor <= 0.5 ) go to 15
     factor = factor*wm1
   10 continue
  ip = 940
   15 ipstr = ip
  return
end
subroutine iptgen (ncolor,ipt,nc)
!
!*******************************************************************************
!
!! IPTGEN generates the pointer vector to block rows.
!
!
!    The algorithm is for block structured matrices with nonconstant block size.
!
!  Parameters:
!
!     ncolor   the number of colors (block rows)
!     ipt      upon input, an integer vector of length ncolor+1
!              upon output, the pointer vector
!     nc       integer vector of length ncolor giving the
!               number of nodes for each color
!
!  
!
  integer ipt(1), nc(1)
!
  ipt(1) = 0
  do 10 k = 1,ncolor
     ipt(k+1) = ipt(k) + nc(k)
 10   continue
  return
end
subroutine itcg (suba,subq,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs,r,p,z,tri,ier)
!
!*******************************************************************************
!
!! ITCG does the conjugate gradient iterations.
!
!
!  Parameters:
!
!         suba      matrix-vector multiplication routine
!         subq      preconditioning routine
!         n         order of system (= nn)
!         u         current solution
!         ubar      known solution (optional)
!         rhs       right hand side vector
!         r,p,z     workspace vectors of length n each
!         tri       tridiagonal matrix associated with the
!                    eigenvalues of the tridiagonal matrix.
!         ier       error code
!
!  
!
  external suba, subq
  integer jcoef(2), jwfac(1)
  dimension coef(1), wfac(1)
  dimension u(1), ubar(1), rhs(1), r(1), p(1), z(1), tri(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
  n = nn
  in = 0
  is = 0
  rzdot = 0.0
  alpha = 0.0
  beta = 0.0
  alphao = 0.0
  maxadp = maxadd
  minadp = minadd
!
!     compute r = residual
!
  call suba (coef,jcoef,wfac,jwfac,n,u,r)
  do 10 i = 1,n
 10   r(i) = rhs(i) - r(i)
  go to 25
!
!************* begin iteration loop.
!
 15   do 20 i = 1,n
 20   r(i) = r(i) - alpha*z(i)
!
!  do preconditioning step -- solve q*z = r for z.
!
 25   call subq (coef,jcoef,wfac,jwfac,n,r,z)
!
!  compute rzdot = (r,z)
!
  dkm1 = rzdot
  rzdot = 0.0
  do 30 i = 1,n
 30   rzdot = rzdot + r(i)*z(i)
  if (rzdot > 0.0) go to 35
  ier = -7
  call ershow (ier,'itcg')
  return
!
!  determine whether or not to stop.
!
 35   call pstops (n,r,z,u,ubar,ier)
  if (level >= 2) call iterm (n,u)
  if (halt  .or.  ier < 0) return
  if (in < itmax) go to 40
  ier = 1
  call ershow (ier,'itcg')
  zeta = stptst
  return
!
!  compute   beta = rzdot/dkm1
!
 40   if (in == 0) go to 45
  beta = rzdot/dkm1
!
!  compute   p = z + beta*p
!
 45   do 50 i = 1,n
 50   p(i) = z(i) + beta*p(i)
!
!  compute   alpha = rzdot / (p,a*p)
!
  call suba (coef,jcoef,wfac,jwfac,n,p,z)
  alphao = alpha
  pap = 0.0
  do 55 i = 1,n
 55   pap = pap + p(i)*z(i)
  alpha = rzdot / pap
  if (pap > 0.0) go to 60
  ier = -6
  call ershow (ier,'itcg')
  return
!
!  compute latest eigenvalue estimates.
!
 60   if (maxadp .or. minadp) call chgcon (tri,ier)
!
!  compute new solution   u = u + alpha*p
!
  do 65 i = 1,n
 65   u(i) = u(i) + alpha*p(i)
  in = in + 1
  is = is + 1
  go to 15
end
subroutine iterm (nn,u)
!
!*******************************************************************************
!
!! ITERM produces the iteration summary line at the end of each iteration. 
!
!
!  if level >= 4, the latest approximation
!     to the solution will be printed.
!
!  Parameters:
!
!          n      order of system  (= nn)
!          u      solution estimate
!
!  
!
  dimension u(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
  n = nn
!
!  print various parameters after each iteration
!
  if (in > 0) go to 15
!
!  print header
!
  if (iacel /= 3) write (nout,10)
 10   format (/5x,'intermediate output after each iteration' &
      /' iteration',11x,'convergence ', &
      5x,'emax',9x,'emin' /7x,'n',7x,'s',8x,'test' /)
  if (iacel == 3) write (nout,12)
 12   format (////5x,'intermediate output after each iteration' &
            //' number of',11x,'convergence',5x, &
            'emax',8x,'omega',7x,'spectral' /' iterations', &
            13x,'test',34x,'radius' //)
!
!  print summary line
!
 15   if (iacel /= 3) write (nout,20) in,is,stptst,emax,emin
 20   format (3x,i5,3x,i5,3x,3e13.5)
  if (iacel == 3) write (nout,22) in,is,stptst,emax,omega,specr
 22   format (3x,i5,3x,i5,3x,5e13.5)
  if (level >= 4) go to 25
  return
!
 25   write (nout,30) in
 30   format (/1x,2x,'estimate of solution at iteration ',i5)
  write (nout,35) (u(i),i=1,n)
 35   format (1x,5g16.7)
  write (nout,40)
 40   format (//)
!
  return
end
subroutine itsi (suba,subq,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs,r,p,z,wksp,ier)
!
!*******************************************************************************
!
!! ITSI does the semi-iterative iterations.
!
!
!  Parameters:
!
!         suba      matrix-vector multiplication routine
!         subq      preconditioning routine
!         n         order of system (= nn)
!         u         current solution
!         ubar      known solution (optional)
!         rhs       right hand side vector
!         r,p,z,    workspace vectors of length n each
!         wksp      volatile workspace
!         ier       error code
!
!  
!
  external suba, subq
  integer jcoef(2), jwfac(1)
  dimension coef(1), wfac(1)
  dimension u(1), ubar(1), rhs(1), r(1), p(1), z(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
  n = nn
  in = 0
!
!  new chebychev sequence.
!
 10   is = 0
  alpha = 0.0
  beta = 0.0
  rho = 1.0
  rzdot = 0.0
  gamma = 2.0/(emax + emin)
  sigma = (emax - emin)/(emax + emin)
  term = sqrt (1.0 - sigma*sigma)
  rr = (1.0 - term)/(1.0 + term)
  maxadp = maxadd
  minadp = minadd
!
!     compute r = residual
!
  call suba (coef,jcoef,wfac,jwfac,n,u,r)
  do 15 i = 1,n
 15   r(i) = rhs(i) - r(i)
  go to 30
!
!************* begin iteration loop.
!
 20   do 25 i = 1,n
 25   r(i) = r(i) - alpha*z(i)
!
!  do preconditioning step -- solve q*z = r for z.
!
 30   call subq (coef,jcoef,wfac,jwfac,n,r,z)
!
!  compute rzdot = (r,z)
!
  dkm1 = rzdot
  rzdot = 0.0
  do 35 i = 1,n
 35   rzdot = rzdot + r(i)*z(i)
  if (is == 0) dkq = rzdot
  if (rzdot >= 0.0) go to 40
  ier = -7
  call ershow (ier,'itsi')
  return
!
!  determine whether or not to stop.
!
 40   call pstops (n,r,z,u,ubar,ier)
  if (level >= 2) call iterm (n,u)
  if (halt  .or.  ier < 0) return
  if (in < itmax) go to 45
  ier = 1
  call ershow (ier,'itsi')
  zeta = stptst
  return
!
!  compute iteration parameters.
!
 45   call parsi
!
!  compute   p = z + beta*p
!            u = u + alpha*p
!
  do 50 i = 1,n
     p(i) = z(i) + beta*p(i)
     u(i) = u(i) + alpha*p(i)
 50   continue
!
!  adapt on emin and emax
!
  in = in + 1
  if (.not. maxadp  .and.  .not. minadp) go to 55
  call chgsi (suba,coef,jcoef,wfac,jwfac,n,z,wksp,icode,ier)
  if (ier < 0) return
!
!  check if new estimates of emax, emin are to be used.
!
  if (icode == 1) go to 10
!
!  estimates of emax, emin are still good.
!
 55   is = is + 1
  call suba (coef,jcoef,wfac,jwfac,n,p,z)
  go to 20
end
subroutine itsor (subq,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs,wksp,ier)
!
!*******************************************************************************
!
!! ITSOR does the SOR iterations.
!
!
!  Parameters:
!
!          subq   routine to do an sor pass
!          n      size of system
!          rhs    right hand side
!          u      solution vector
!          ubar   known solution (optional)
!          wksp   workspace vector of length 2*n
!
!  
!
  integer jcoef(2), jwfac(1)
  dimension coef(1), wfac(1)
  dimension rhs(1), u(1), ubar(1), wksp(1)
  external subq
  logical  change
!
! *** begin -- itpack common
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
!
!  set initial parameters not already set
!
  n = nn
  in = 0
  is = 0
  ip = 0
  iss = 0
  iphat = 2
  delnnm = 0.0
  delsnm = 0.0
  call sorstp (n,u,ubar,0.0,0.0)
  change = omgadp
  ib2 = n + 1
  if (.not. omgadp) go to 10
  omegap = omega
  omega = 1.0
  ipstar = 4
  if (omegap <= 1.0) change = .false.
!
!  start iterating.
!
 10   do 55 iter = 1,itmax+1
!
!  output intermediate information
!
     if (level >= 2) call iterm (n,u)
     if (halt) return
     if (.not. change) go to 15
     change = .false.
     is = is + 1
     ip = 0
     iss = 0
     omega = amin1 (omegap,tau(is))
     iphat = max ( 3 , int ( (omega-1.0)/(2.0-omega) ) )
     ipstar = ipstr (omega)
!
!  compute u (in + 1) and norm of del(s,p)
!
 15      delsnm = delnnm
     spcrm1 = specr
     do 20 i = 1,n
 20      wksp(i) = rhs(i)
     call subq (coef,jcoef,wfac,jwfac,n,u,wksp,wksp(ib2))
     do 25 i = 1,n
 25      wksp(i) = u(i) - wksp(n+i)
     sum = 0.0
     do 28 i = 1,n
 28      sum = sum + wksp(i)*wksp(i)
     delnnm = sqrt (sum)
     do 30 i = 1,n
 30      u(i) = wksp(i+n)
     if (delnnm == 0.0) go to 35
     if (in /= 0) specr = delnnm / delsnm
     if (ip < iphat) go to 50
!
!  stopping test, set h
!
     if (specr >= 1.0) go to 50
     if (.not. (specr > (omega - 1.0))) go to 35
     h = specr
     go to 40
 35      iss = iss + 1
     h = omega - 1.0
!
!  perform stopping test.
!
 40      continue
     dnrm = delnnm**2
     call sorstp (n,u,ubar,dnrm,h)
     if (halt) go to 50
!
!  method has not converged yet, test for changing omega
!
     if (.not. omgadp) go to 50
     if (ip < ipstar)  go to 50
     if (omega > 1.0) go to 45
     emax = sqrt (abs (specr))
     omegap = 2.0 / (1.0 + sqrt (abs (1.0 - specr)))
     change = .true.
     go to 50
 45      if (iss /= 0) go to 50
     if (specr <= (omega - 1.0)**fff) go to 50
     if ((specr + 0.00005) <= spcrm1) go to 50
!
!  change parameters
!
     emax = (specr + omega - 1.0) / (sqrt (abs (specr))*omega)
     omegap = 2.0 / (1.0 + sqrt (abs (1.0 - emax*emax)))
     change = .true.
!
 50      ip = ip + 1
     in = in + 1
 55   continue
  ier = 1
  in = in - 1
  call ershow (ier,'itsor')
  zeta = stptst
  return
end
subroutine itsrcg (suba,subq,subadp,coef,jcoef,wfac,jwfac,nn,u,ubar, &
  rhs,r,p,z,tri,ier)
!
!*******************************************************************************
!
!! ITSRCG does the SSOR conjugate gradient iterations.
!
!
!  Parameters:
!
!         suba      matrix-vector multiplication routine
!         subq      preconditioning routine
!         subadp    adpation routine
!         n         order of system (= nn)
!         u         current solution
!         ubar      known solution (optional)
!         rhs       right hand side vector
!         r,p,z     workspace vectors of length n each
!         tri       tridiagonal matrix associated with the
!                    eigenvalues of the tridiagonal matrix.
!         ier       error code
!
!  
!
  external suba, subq, subadp
  integer jcoef(2), jwfac(1)
  dimension coef(1), wfac(1)
  dimension u(1), ubar(1), rhs(1), r(1), p(1), z(1), tri(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
  n = nn
  in = 0
  isw = 1
 5    is = 0
  rzdot = 0.0
  alpha = 0.0
  beta = 0.0
  alphao = 0.0
  maxadp = maxadd
  minadp = minadd
!
!     recompute bnorm
!
  call nmcalc (coef,jcoef,wfac,jwfac,isw,subq,n,rhs,ubar,r,ier)
  if (ier < 0) return
  isw = 2
!
!     compute r = residual
!
  call suba (coef,jcoef,wfac,jwfac,n,u,r)
  do 10 i = 1,n
 10   r(i) = rhs(i) - r(i)
  go to 25
!
!************* begin iteration loop.
!
 15   do 20 i = 1,n
 20   r(i) = r(i) - alpha*z(i)
!
!  do preconditioning step -- solve q*z = r for z.
!
 25   call subq (coef,jcoef,wfac,jwfac,n,r,z)
!
!  compute rzdot = (r,z)
!
  dkm1 = rzdot
  rzdot = 0.0
  do 30 i = 1,n
 30   rzdot = rzdot + r(i)*z(i)
  if (rzdot >= 0.0) go to 35
  ier = -7
  call ershow (ier,'itsrcg')
  return
!
!  determine whether or not to stop.
!
 35   call pstops (n,r,z,u,ubar,ier)
  if (level >= 2) call iterm (n,u)
  if (halt  .or.  ier < 0) return
  if (in < itmax) go to 40
  ier = 1
  call ershow (ier,'itsrcg')
  zeta = stptst
  return
!
!  compute   beta = rzdot/dkm1
!
 40   if (is == 0) go to 45
  beta = rzdot/dkm1
!
!  compute   p = z + beta*p
!
 45   do 50 i = 1,n
 50   p(i) = z(i) + beta*p(i)
!
!  compute   alpha = rzdot / (p,a*p)
!
  call suba (coef,jcoef,wfac,jwfac,n,p,z)
  alphao = alpha
  pap = 0.0
  do 55 i = 1,n
 55   pap = pap + p(i)*z(i)
  alpha = rzdot / pap
  if (pap >= 0.0) go to 60
  ier = -6
  call ershow (ier,'itsrcg')
  return
!
!  compute latest eigenvalue estimates.
!
 60   if (minadp) call chgcon (tri,ier)
!
!  compute new solution   u = u + alpha*p
!
  do 65 i = 1,n
 65   u(i) = u(i) + alpha*p(i)
  is = is + 1
  in = in + 1
  call ssorad (subadp,coef,jcoef,wfac,jwfac,n,p,z,r,icode)
  if (icode == 0) go to 15
  go to 5
end
subroutine itsrsi (suba,subq,subadp,coef,jcoef,wfac,jwfac,nn,u,ubar, &
  rhs,r,p,z,wksp,ier)
!
!*******************************************************************************
!
!! ITSRSI does the SSOR semi-iterative iterations.
!
!  Parameters:
!
!         suba      matrix-vector multiplication routine
!         subq      preconditioning routine
!         subadp    adpation routine
!         n         order of system (= nn)
!         u         current solution
!         ubar      known solution (optional)
!         rhs       right hand side vector
!         r,p,z,    workspace vectors of length n each
!         wksp      volatile workspace
!         ier       error code
!
!  
!
  external suba, subq, subadp
  integer jcoef(2), jwfac(1)
  dimension coef(1), wfac(1)
  dimension u(1), ubar(1), rhs(1), r(1), p(1), z(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
  n = nn
!
  in = 0
  isw = 1
!
!     recompute bnorm
!
 5    call nmcalc (coef,jcoef,wfac,jwfac,isw,subq,n,rhs,ubar,r,ier)
  if (ier < 0) return
  isw = 2
!
!  update rayleigh quotient .
!
  if (in == 0) go to 10
  call subq (coef,jcoef,wfac,jwfac,n,p,z)
  call suba (coef,jcoef,wfac,jwfac,n,z,r)
  rq = vdot (n,z,r) / vdot (n,z,p)
  rqmin = rq
  if (minadd) emin = rqmin
!
!  new chebychev sequence.
!
 10   is = 0
  alpha = 0.0
  beta = 0.0
  rho = 1.0
  rzdot = 0.0
  gamma = 2.0/(emax + emin)
  sigma = (emax - emin)/(emax + emin)
  term = sqrt (1.0 - sigma*sigma)
  rr = (1.0 - term)/(1.0 + term)
  minadp = minadd
!
!     compute r = residual
!
  call suba (coef,jcoef,wfac,jwfac,n,u,r)
  do 15 i = 1,n
 15   r(i) = rhs(i) - r(i)
  go to 30
!
!************* begin iteration loop.
!
 20   do 25 i = 1,n
 25   r(i) = r(i) - alpha*z(i)
!
!  do preconditioning step -- solve q*z = r for z.
!
 30   call subq (coef,jcoef,wfac,jwfac,n,r,z)
!
!  compute rzdot = (r,z)
!
  dkm1 = rzdot
  rzdot = 0.0
  do 35 i = 1,n
 35   rzdot = rzdot + r(i)*z(i)
  if (is == 0) dkq = rzdot
  if (rzdot >= 0.0) go to 40
  ier = -7
  call ershow (ier,'itsrsi')
  return
!
!  determine whether or not to stop.
!
 40   call pstops (n,r,z,u,ubar,ier)
  if (level >= 2) call iterm (n,u)
  if (halt  .or.  ier < 0) return
  if (in < itmax) go to 45
  ier = 1
  call ershow (ier,'itsrsi')
  zeta = stptst
  return
!
!  compute iteration parameters.
!
 45   call parsi
!
!  compute   p = z + beta*p
!            u = u + alpha*p
!
  do 50 i = 1,n
     p(i) = z(i) + beta*p(i)
     u(i) = u(i) + alpha*p(i)
 50   continue
!
!  adapt on emin and emax
!
  in = in + 1
  if (.not. minadp) go to 55
  call chgsi (suba,coef,jcoef,wfac,jwfac,n,z,wksp,icode,ier)
  if (ier < 0) return
!
!  check if new estimates of emax, emin are to be used.
!
  if (icode == 1) go to 10
!
!  estimates of emax, emin are still good.
!
 55   is = is + 1
  call suba (coef,jcoef,wfac,jwfac,n,p,z)
  call ssorad (subadp,coef,jcoef,wfac,jwfac,n,p,z,r,icode)
  if (icode == 0) go to 20
  go to 5
end
subroutine jac1 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! JAC1 drives the Jacobi preconditioner.
!
  external accel, suba8, suba9, subq1, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom4 / keygs, srelpr, keyzer
!
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + n
  call split (accel,suba8,suba9,subq1,subq1,subq1,subq1,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  if (keygs == 1) irpnt = irpnt - n
  return
end
subroutine jac2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! JAC2 drives the Jacobi preconditioner.
!
  external accel, suba1, subq1, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  call split (accel,suba1,suba1,subq1,subq1,subq1,subq1,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine jac3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! JAC3 drives the Jacobi preconditioner.
!
  external accel, suba4, suba5, subq1, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  call split (accel,suba4,suba5,subq1,subq1,subq1,subq1,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine jac4 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! JAC4 drives the Jacobi preconditioner.
!
  external accel, suba12, subq1, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom4 / keygs, srelpr, keyzer
!
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba12,suba12,subq1,subq1,subq1,subq1,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine jac5 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! JAC5 drives the Jacobi preconditioner.
!
  external accel, suba13, suba14, subq1, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom4 / keygs, srelpr, keyzer
!
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba13,suba14,subq1,subq1,subq1,subq1,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine landir (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef, &
  n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LANDIR is the user interface to the Lanczos/ORTHODIR algorithm. 
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call ldirw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wksp,iwksp,n, &
    u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine lanmin (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n, &
  u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LANMIN is the user interface to the Lanczos/ORTHOMIN algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call lminw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wksp,iwksp,n, &
    u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine lanres (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n, &
  u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LANRES is the user interface to the Lanczos/ORTHORES algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call lresw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wksp,iwksp,n, &
    u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine ldirw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LDIRW runs the Lanczos/ORTHODIR algorithm. 
!
!  see jea and young, in
! linear algebra and its applications, vol 52/3, 1983, pp399f.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subat, subql, subqlt, subqr, subqrt
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! indexing functions.
!
  indq(i) = iq + n*mod(i,2)
  indqt(i) = iqt + n*mod(i,2)
!
! preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 14
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  if (iqr) go to 995
  if (level >= 2) write (nout,496)
496   format (' landir')
!
! initialize the stopping test.
!
  call inithv (0)
  zhave  = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
  iq = 1
  iqt = iq + 2*n
  ir = iqt + 2*n
  iv1 = ir + n
  iv2 = iv1 + n
  iv3 = iv2 + n
  nwusd = max(nwusd,iv3-1+n)
!
! check the memory usage.
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(ir))
!
! begin iteration loop
!
! determine whether or not to stop.
!
 10   call inithv (1)
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,wk(ir),xxx,wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
  if (in /= 0) go to 110
!
! perform first-iterate calculations
!
  call vcopy (n,wk(ir),wk(indq(in)))
  call vcopy (n,wk(indq(in)),wk(indqt(in)))
  qaq= 0.0
  go to 115
!
! proceed to calculate the direction vectors, for in > 0.
!
 110  call subqlt (coef,jcoef,wfac,jwfac,n,wk(indqt(in-1)),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv3))
  aqaq = vdot(n,wk(iv2),wk(iv3))
  an = aqaq / qaq
  if (in /= 1) go to 150
  call vtriad (n,wk(indq(in)),wk(iv2),-an,wk(indq(in-1)),1)
  call vtriad (n,wk(indqt(in)),wk(iv3),-an,wk(indqt(in-1)),1)
  go to 115
 150  bn = qaq / qaqold
  call vtriad (n,wk(indq(in)),wk(iv2),-bn,wk(indq(in-2)),1)
  call vtriad (n,wk(indq(in)),wk(indq(in)),-an,wk(indq(in-1)),1)
  call vtriad (n,wk(indqt(in)),wk(iv3),-bn,wk(indqt(in-2)),1)
  call vtriad (n,wk(indqt(in)),wk(indqt(in)),-an,wk(indqt(in-1)),1)
!
! proceed to form the iterate.
!
 115  call suba (coef,jcoef,wfac,jwfac,n,wk(indq(in)),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  qaqold = qaq
  qaq = vdot(n,wk(indqt(in)),wk(iv2))
  if (abs(qaq) < srelpr) go to 998
  qr = vdot(n,wk(indqt(in)),wk(ir))
  vlamda = qr / qaq
  call vtriad (n,u,u,vlamda,wk(indq(in)),1)
  call vtriad (n,wk(ir),wk(ir),-vlamda,wk(iv2),1)
!
! proceed to next iteration
!
  in = in + 1
  is = is + 1
  go to 10
!
!  Finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'ldirw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' lanczos/orthodir converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 995  ier = -16
  call ershow (ier,'ldirw')
  return
!
 997  call ershow (ier,'ldirw')
  go to 735
!
 998  ier = -15
  call ershow (ier,'ldirw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'ldirw')
  go to 735
!
end
subroutine lfact (coef,jcoef,wksp,nn,ier)
!
!*******************************************************************************
!
!! LFACT computes a line factorization.
!
!
!  Parameters:
!
!        n        problem size
!        nfactr   factorization size
!
!  common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2)
  dimension coef(1), wksp(1)
!
!  check for sufficient workspace to store factor.
!
  n = nn
  if (nstore == 2) isym = 0
  if (nstore == 3) isym = 1
  ndt = 0
  ndb = 0
  do 20 jd = 1,maxnz
     do 15 j = 1,maxnz
        if (jcoef(j) /= jd) go to 15
        ndt = ndt + 1
        go to 20
 15      continue
     go to 25
 20   continue
 25   if (isym == 0) go to 40
  do 35 jd = 1,maxnz
     do 30 j = 1,maxnz
        if (jcoef(j) /= -jd) go to 30
        ndb = ndb + 1
        go to 35
 30      continue
     go to 40
 35   continue
 40   nfactr = (ndt + ndb + 1)*n
  call needw ('lfact',0,irpnt,nfactr,ier)
  if (ier < 0) return
!
  ifactr = irpnt
  call vcopy (n,coef,wksp(ifactr))
  ndt = 0
  do 55 jd = 1,maxnz
     do 50 j = 1,maxnz
        if (jcoef(j) /= jd) go to 50
        ndt = ndt + 1
        ipt1 = (j - 1)*ndim + 1
        ipt2 = ndt*n + ifactr
        call vcopy (n,coef(ipt1),wksp(ipt2))
        go to 55
 50      continue
     go to 60
 55   continue
 60   ndb = 0
  if (isym == 0) go to 75
  do 70 jd = 1,maxnz
     do 65 j = 1,maxnz
        if (jcoef(j) /= -jd) go to 65
        ndb = ndb + 1
        ipt1 = (j - 1)*ndim + 1
        ipt2 = (ndt + ndb)*n + ifactr
        call vcopy (n,coef(ipt1),wksp(ipt2))
        go to 70
 65      continue
     go to 75
 70   continue
!
!  factor.
!
 75   call bdfac (n,n,kblsz,ndt,ndb,wksp(ifactr),isym)
  irpnt = irpnt + nfactr
  return
end
subroutine linv (coef,jcoef,wksp,nn,ier)
!
!*******************************************************************************
!
!! LINV computes a line approximate inverse.
!
!
!  Parameters:
!
!        n        problem size
!        nfactr   factorization size
!
!  common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2)
  dimension coef(1), wksp(1)
!
!  check for sufficient workspace to store factor.
!
  n = nn
  if (nstore == 2) isym = 0
  if (nstore == 3) isym = 1
  ndt = 0
  ndb = 0
  do 20 jd = 1,maxnz
     do 15 j = 1,maxnz
        if (jcoef(j) /= jd) go to 15
        ndt = ndt + 1
        go to 20
 15      continue
     go to 25
 20   continue
 25   if (isym == 0) go to 40
  do 35 jd = 1,maxnz
     do 30 j = 1,maxnz
        if (jcoef(j) /= -jd) go to 30
        ndb = ndb + 1
        go to 35
 30      continue
     go to 40
 35   continue
!
 40   ndt = ndt + ltrunc
  if (isym == 1) ndb = ndb + ltrunc
  nfactr = (ndt + ndb + 1)*n
  call needw ('linv',0,irpnt,nfactr,ier)
  if (ier < 0) return
!
  ifactr = irpnt
  call vfill (nfactr,wksp(ifactr),0.0)
  call vcopy (n,coef,wksp(ifactr))
  it = 0
  do 55 jd = 1,maxnz
     do 50 j = 1,maxnz
        if (jcoef(j) /= jd) go to 50
        it = it + 1
        ipt1 = (j - 1)*ndim + 1
        ipt2 = it*n + ifactr
        call vcopy (n,coef(ipt1),wksp(ipt2))
        go to 55
 50      continue
     go to 60
 55   continue
 60   if (isym == 0) go to 75
  it = ndt
  do 70 jd = 1,maxnz
     do 65 j = 1,maxnz
        if (jcoef(j) /= -jd) go to 65
        it = it + 1
        ipt1 = (j - 1)*ndim + 1
        ipt2 = it*n + ifactr
        call vcopy (n,coef(ipt1),wksp(ipt2))
        go to 70
 65      continue
     go to 75
 70   continue
!
!  factor and invert.
!
 75   call bdfac (n,n,kblsz,ndt,ndb,wksp(ifactr),isym)
  call bdinv (n,n,kblsz,ndt,ndb,wksp(ifactr),isym)
  irpnt = irpnt + nfactr
  return
end
subroutine ljac2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LJAC2 drives the line Jacobi preconditioner.
!
  external accel, suba1, subq2, noadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call split (accel,suba1,suba1,subq2,subq2,subq2,subq2,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine ljac3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LJAC3 drives the line Jacobi preconditioner.
!
  external accel, suba4, suba5, subq2, subq3, noadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call split (accel,suba4,suba5,subq2,subq3,subq2,subq3,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine ljacx2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LJACX2 drives the line Jacobi preconditioner.
!
  external accel, suba1, subq4, noadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  t1 = timer (dummy)
  if (ifact == 1) call linv (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call split (accel,suba1,suba1,subq4,subq4,subq4,subq4,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine ljacx3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LJACX3 drives the line Jacobi preconditioner.
!
  external accel, suba4, suba5, subq4, subq5, noadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  t1 = timer (dummy)
  if (ifact == 1) call linv (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call split (accel,suba4,suba5,subq4,subq5,subq4,subq5,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine llsp2 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LLSP2 drives the line least squares polynomial preconditioner.
!
  external accel, suba1, subq23, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom8 / ainf
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  n = nn
  call needw ('llsp2',0,irpnt,n,ier)
  if (ier < 0) return
  call adinfn (n,ndim,maxnz,jcoef,coef,2,ainf,wksp(irpnt))
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call needw ('llsp2',0,irpnt,2*n,ier)
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*n
  call split (accel,suba1,suba1,subq23,subq23,subq23,subq23,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  return
end
subroutine llsp3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LLSP3 drives the line least squares polynomial preconditioner.
!
  external accel, suba4, suba5, subq66, subq67, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom8 / ainf
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  n = nn
  call needw ('llsp3',0,irpnt,n,ier)
  if (ier < 0) return
  call adinfn (n,ndim,maxnz,jcoef,coef,3,ainf,wksp(irpnt))
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call needw ('llsp3',0,irpnt,2*n,ier)
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*n
  call split (accel,suba4,suba5,subq66,subq67,subq66,subq67,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  return
end
subroutine lminw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LMINW runs the Lanczos/ORTHOMIN algorithm.
!
!
! here, zhat and phat will refer to the "dummy" system of the
! lanczos method.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2),wfac(1), jwfac(1)
  integer vect1, vect2, os
  external suba, subat, subql, subqlt, subqr, subqrt
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
!
  nwusd = 0
  ier = 0
  iacel = 15
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  if (level >= 2) write (nout,496)
496   format (' lanmin')
!
! initialize the stopping test.
!
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  call inithv (0)
  zhave  = .true.
  zthave = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
! allocate memory -- overlap wherever possible.
  ip = 1
  ipt = ip + n
  if (.not. iqr) ipt = ip
  iphat = ipt + n
  iz = iphat + n
  izt = iz + n
  if (.not. iqr) izt = iz
  izhat = izt + n
  iv1 = izhat + n
  iv2 = iv1 + n
  if (iqlr == 0) nwusd = max(nwusd,iv1-1+n)
  if (iqlr /= 0) nwusd = max(nwusd,iv2-1+n)
!
! check the memory usage.
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  if (.not. iql) go to 121
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iz))
  go to 122
 121  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iz))
  call vexopy (n,wk(iz),rhs,wk(iz),2)
 122  if (iqr) call subqr (coef,jcoef,wfac,jwfac,n,wk(iz),wk(izt))
!
!  begin iteration loop
!
! determine whether or not to stop.
!
 10   call inithv (1)
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,wk(iz),wk(izt),wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
  if (in /= 0) go to 110
!
! perform first-iterate calculations
!
  call vcopy (n,wk(iz),wk(izhat))
  rd = vdot (n,wk(iz),wk(izhat))
  call vcopy (n,wk(iz),wk(ip))
  call vcopy (n,wk(izhat),wk(iphat))
  if (iqr) call vcopy (n,wk(izt),wk(ipt))
  go to 111
!
! perform subsequent-iterate calculations
!
 110  rdold = rd
!     if (abs(rdold) < srelpr) go to 996
  if (abs(rdold) == 0e0) go to 996
!
! form the old zhat.
  go to (131,132,133,134), iqlr + 1
 131  call subat (coef,jcoef,wfac,jwfac,n,wk(iphat),wk(iv1))
  go to 135
 132  call subqlt (coef,jcoef,wfac,jwfac,n,wk(iphat),wk(iv2))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv2),wk(iv1))
  go to 135
 133  call subat (coef,jcoef,wfac,jwfac,n,wk(iphat),wk(iv2))
  call subqrt (coef,jcoef,wfac,jwfac,n,wk(iv2),wk(iv1))
  go to 135
 134  call subqlt (coef,jcoef,wfac,jwfac,n,wk(iphat),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  call subqrt (coef,jcoef,wfac,jwfac,n,wk(iv2),wk(iv1))
 135  call vtriad (n,wk(izhat),wk(izhat),-vlamda,wk(iv1),1)
!
! form the direction vectors.
  rd = vdot (n,wk(iz),wk(izhat))
  an = rd/rdold
  call vtriad (n,wk(ip),wk(iz),an,wk(ip),1)
  call vtriad (n,wk(iphat),wk(izhat),an,wk(iphat),1)
  if (iqr) call vtriad (n,wk(ipt),wk(izt),an,wk(ipt),1)
!
!  Form the iterate.
!
 111  if (iql) go to 141
  call suba (coef,jcoef,wfac,jwfac,n,wk(ipt),wk(iv1))
  go to 142
 141  call suba (coef,jcoef,wfac,jwfac,n,wk(ipt),wk(iv2))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv2),wk(iv1))
 142  pap = vdot (n,wk(iphat),wk(iv1))
!     if (abs(pap) < srelpr**2) go to 998
  if (abs(pap) == 0e0) go to 998
  vlamda = rd/pap
!
!
  call vtriad (n,u,u,vlamda,wk(ipt),1)
  call vtriad (n,wk(iz),wk(iz),-vlamda,wk(iv1),1)
  if (.not. iqr) go to 151
  call subqr (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  call vtriad (n,wk(izt),wk(izt),-vlamda,wk(iv2),1)
!
! proceed to next iteration
!
 151  in = in + 1
  is = is + 1
  go to 10
!
!  Finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'lminw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' lanczos/orthomin converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 996  ier = -13
  call ershow (ier,'lminw')
  go to 725
!
 997  call ershow (ier,'lminw')
  go to 735
!
 998  ier = -15
  call ershow (ier,'lminw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'lminw')
  go to 735
!
end
subroutine lneu2 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LNEU2 drives the line Neumann polynomial preconditioner.
!
  external accel, suba1, subq24, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  n = nn
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call needw ('lneu2',0,irpnt,2*n,ier)
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*n
  call split (accel,suba1,suba1,subq24,subq24,subq24,subq24,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  return
end
subroutine lneu3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LNEU3 drives the line Neumann polynomial preconditioner.
!
  external accel, suba4, suba5, subq68, subq69, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  n = nn
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call needw ('lneu3',0,irpnt,2*n,ier)
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*n
  call split (accel,suba4,suba5,subq68,subq69,subq68,subq69,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  return
end
subroutine lresw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LRESW runs the Lanczos/ORTHORES algorithm.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subat, subql, subqlt, subqr, subqrt
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! indexing functions.
!
  indu(i) = iu + n*mod(i,nv)
  indr(i) = ir + n*mod(i,nv)
  indrt(i) = irt + n*mod(i,nv)
!
! preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 16
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  if (iqr) go to 995
  if (level >= 2) write (nout,496)
496   format (' lanres')
!
! initialize the stopping test.
!
  call inithv (0)
  zhave  = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
  nv = 2
  iu = 1
  ir = iu + nv*n
  irt = ir + nv*n
  iv1 = irt + nv*n
  nwusd = max(nwusd,iv1-1+n)
!
! check the memory usage.
!
  if (nwusd > nw) go to 999
!
! note -- we will use the vector 'u' for scratch storage, to save space.
!
  call vcopy (n,u,wk(indu(0)))
  in = 0
  is = 0
  call suba (coef,jcoef,wfac,jwfac,n,wk(indu(in)),wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(indr(in)))
  call vcopy (n,wk(indr(in)),wk(indrt(in)))
!
!  Begin iteration loop.
!
! determine whether or not to stop.
!
 10   call inithv (1)
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,wk(indu(in)), &
    ubar,rhs,xxx,wk(indr(in)),xxx,wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,wk(indu(in)))
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
! proceed to calculate the parameters.
! first, gamma.
!
  rd = vdot (n,wk(indr(in)),wk(indrt(in)))
  call suba (coef,jcoef,wfac,jwfac,n,wk(indr(in)),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),u)
  rar = vdot (n,u,wk(indrt(in)))
  if (abs(rar) < srelpr) go to 998
  gam = rd / rar
!
! now, rho.
!
  if (in /= 0) go to 118
  rho = 1.0
  go to 119
 118  if (abs(gamold) < srelpr) go to 998
  if (abs(rdold) < srelpr) go to 998
  if (abs(rho) < srelpr) go to 998
  rhoinv = 1.0 - (gam/gamold)*(rd/rdold)/rho
  if (abs(rhoinv) < srelpr) go to 998
  rho = 1.0 / rhoinv
!
! now work on updating u, r, rt.
! first, the first iteration.
!
 119  if (in /= 0) go to 150
  call vtriad (n,wk(indu(in+1)),wk(indu(in)),gam,wk(indr(in)),1)
  call vtriad (n,wk(indu(in+1)),xxx,rho,wk(indu(in+1)),2)
  call vtriad (n,wk(indr(in+1)),wk(indr(in)),-gam,u,1)
  call vtriad (n,wk(indr(in+1)),xxx,rho,wk(indr(in+1)),2)
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(indrt(in)),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),u)
  call vtriad (n,wk(indrt(in+1)),wk(indrt(in)),-gam,u,1)
  call vtriad (n,wk(indrt(in+1)),xxx,rho,wk(indrt(in+1)),2)
  go to 151
!
! now work on subsequent iterations.
!
 150  call vtriad (n,wk(indu(in+1)),xxx,1.0-rho,wk(indu(in-1)),2)
  call vtriad (n,wk(indu(in+1)),wk(indu(in+1)),rho,wk(indu(in)),1)
  call vtriad (n,wk(indu(in+1)),wk(indu(in+1)),rho*gam,wk(indr(in)),1)
  call vtriad (n,wk(indr(in+1)),xxx,1.0-rho,wk(indr(in-1)),2)
  call vtriad (n,wk(indr(in+1)),wk(indr(in+1)),rho,wk(indr(in)),1)
  call vtriad (n,wk(indr(in+1)),wk(indr(in+1)),-rho*gam,u,1)
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(indrt(in)),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),u)
  call vtriad (n,wk(indrt(in+1)),xxx,1.0-rho,wk(indrt(in-1)),2)
  call vtriad (n,wk(indrt(in+1)),wk(indrt(in+1)),rho,wk(indrt(in)),1)
  call vtriad (n,wk(indrt(in+1)),wk(indrt(in+1)),-rho*gam,u,1)
!
! proceed to next iteration
!
 151  gamold = gam
  rdold = rd
  in = in + 1
  is = is + 1
  go to 10
!
!  Finish up.
!
 900  call vcopy (n,wk(indu(in)),u)
  if (halt) go to 715
  ier = 1
  call ershow (ier,'lresw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' lanczos/orthores converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 995  ier = -16
  call ershow (ier,'lresw')
  return
!
 997  call ershow (ier,'lresw')
  go to 735
!
 998  ier = -15
  call ershow (ier,'lresw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'lresw')
  go to 735
!
end
subroutine lsor2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSOR2 drives the line SOR method.
!
  external accel, suba1, subq20, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call split (accel,suba1,suba1,subq20,subq20,subq20,subq20,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine lsor3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSOR3 drives the line SOR method.
!
  external accel, suba4, suba5, subq58, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call split (accel,suba4,suba5,subq58,subq58,subq58,subq58,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine lsp1 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSP1 drives the least squares polynomial preconditioner.
!
  external accel, suba8, suba9, subq92, subq93, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom8 / ainf
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom4 / keygs, srelpr, keyzer
!
  n = nn
  call needw ('lsp1',0,irpnt,2*n,ier)
  if (ier < 0) return
  call ainfn (n,ndim,maxnz,jcoef,coef,1,ainf,wksp(irpnt))
  iwkpt2 = irpnt
  irpnt = irpnt + 2*n
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + n
  call split (accel,suba8,suba9,subq92,subq93,subq92,subq93,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  if (keygs == 1) irpnt = irpnt - n
  return
end
subroutine lsp2 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSP2 drives the least squares polynomial preconditioner.
!
  external accel, suba1, subq18, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom8 / ainf
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  call needw ('lsp2',0,irpnt,2*n,ier)
  if (ier < 0) return
  call ainfn (n,ndim,maxnz,jcoef,coef,2,ainf,wksp(irpnt))
  iwkpt1 = irpnt
  irpnt = irpnt + 2*n
  call split (accel,suba1,suba1,subq18,subq18,subq18,subq18,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  return
end
subroutine lsp3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSP3 drives the least squares polynomial preconditioner.
!
  external accel, suba4, suba5, subq54, subq55, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom8 / ainf
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  call needw ('lsp3',0,irpnt,2*n,ier)
  if (ier < 0) return
  call ainfn (n,ndim,maxnz,jcoef,coef,3,ainf,wksp(irpnt))
  iwkpt1 = irpnt
  irpnt = irpnt + 2*n
  call split (accel,suba4,suba5,subq54,subq55,subq54,subq55,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  return
end
subroutine lsp4 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSP4 drives the least squares polynomial preconditioner.
!
  external accel, suba12, sub110, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom8 / ainf
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom4 / keygs, srelpr, keyzer
!
  n = nn
  call needw ('lsp4',0,irpnt,2*n,ier)
  if (ier < 0) return
  call ainfn (n,ndim,maxnz,jcoef,coef,4,ainf,wksp(irpnt))
  iwkpt2 = irpnt
  irpnt = irpnt + 2*n
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba12,suba12,sub110,sub110,sub110,sub110,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine lsp5 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSP5 drives the least squares polynomial preconditioner.
!
  external accel, suba13, suba14, sub112, sub113, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom8 / ainf
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom4 / keygs, srelpr, keyzer
!
  n = nn
  call needw ('lsp5',0,irpnt,2*n,ier)
  if (ier < 0) return
  call ainfn (n,ndim,maxnz,jcoef,coef,5,ainf,wksp(irpnt))
  iwkpt2 = irpnt
  irpnt = irpnt + 2*n
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba13,suba14,sub112,sub113,sub112,sub113,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*n
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine lsqr (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSQR is the user interface to the LSQR algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call lsqrw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wksp,iwksp,n,u, &
    ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine lsqrw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSQRW runs the LSQR algorithm.  
!
!
!  the algorithm is taken from
! the article 'lsqr -- an algorithm for sparse linear equations
! and sparse least squares.'
! by c. c. paige amd m. a. saunders, in acm transactions on
! mathematical software, vol. 8, no. 1, march 1982, pp. 43-71.
! the iterates produced are the same as those of cgnr, in exact
! arithmetic, but this should be more stable.  only left
! preconditioning is currently implemented.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  integer vect1, vect2, os
  external suba, subat, subql, subqlt, subqr, subqrt
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 6
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 996
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  if (iqr) go to 995
  if (level >= 2) write (nout,496)
496   format (' lsqr')
!
! initialize the stopping test.
!
  call inithv (0)
  zdhav = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 735
!
!  associated integer variables.
!
  iu = 1
  iv = iu + n
  iw = iv + n
  iv1 = iw + n
  iv2 = iv1 + n
  nwusd = max(nwusd,iv2-1+n)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
!
! now, perform first-iterate calculations
!
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  beta = sqrt(vdot (n,wk(iv2),wk(iv2)))
  if (abs(beta) < srelpr) go to 997
  call vtriad (n,wk(iu),xxx,1.0/beta,wk(iv2),2)
!
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(iu),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  alpha = sqrt(vdot (n,wk(iv2),wk(iv2)))
  if (abs(alpha) < srelpr) go to 997
  call vtriad (n,wk(iv),xxx,1.0/alpha,wk(iv2),2)
!
  call vcopy (n,wk(iv),wk(iw))
  phibar = beta
  rhobar = alpha
  zdot = phibar**2
! if u(0) is zero, then the norm of u(n) can be calculated for free.
! otherwise, i don't know.
!
!  Begin iteration loop.
!
! determine whether or not to stop --
!
 10   call inithv (1)
  zdhav = .true.
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
!  compute the lanczos vectors.
!
  call suba (coef,jcoef,wfac,jwfac,n,wk(iv),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  call vtriad (n,wk(iu),wk(iv2),-alpha,wk(iu),1)
  beta = sqrt(vdot (n,wk(iu),wk(iu)))
  if (abs(beta) < srelpr) go to 997
  call vtriad (n,wk(iu),xxx,1.0/beta,wk(iu),2)
!
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(iu),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  call vtriad (n,wk(iv),wk(iv2),-beta,wk(iv),1)
  alpha = sqrt(vdot (n,wk(iv),wk(iv)))
  if (abs(alpha) < srelpr) go to 997
  call vtriad (n,wk(iv),xxx,1.0/alpha,wk(iv),2)
!
! continue by calculating various scalars.
!
  rho = sqrt(rhobar**2+beta**2)
  if (rho < srelpr) go to 998
  c = rhobar/rho
  s = beta/rho
  theta = s*alpha
  rhobar = -c*alpha
  phi = c*phibar
  phibar = s*phibar
!
! now generate the new u and w vectors.
!
  call vtriad (n,u,u,phi/rho,wk(iw),1)
  call vtriad (n,wk(iw),wk(iv),-theta/rho,wk(iw),1)
!
! proceed to next iteration
!
  zdot = phibar**2
  in = in + 1
  is = is + 1
  go to 10
!
!  Finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'lsqrw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' lsqr converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 995  ier = -16
  call ershow (ier,'lsqrw')
  return
!
 996  call ershow (ier,'lsqrw')
  go to 735
!
 997  ier = -13
  call ershow (ier,'lsqrw')
  go to 725
!
 998  ier = -14
  call ershow (ier,'lsqrw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'lsqrw')
  go to 735
!
end
subroutine lssor2 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSSOR2 drives the line SSOR method.
!
  external accel, suba1, subq21, subq22, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  n = nn
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  iwkpt1 = irpnt
  irpnt = irpnt + n
  if (ier < 0) return
  call split (accel,suba1,suba1,subq21,subq21,subq21,subq21,copy,copy,subq22, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine lssor3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! LSSOR3 drives the line SSOR method.
!
  external accel, suba4, suba5, subq59, subq60, subq61, subq62, subq63
  external subq64, subq65
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  n = nn
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call lfact (coef,jcoef,wksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  iwkpt1 = irpnt
  irpnt = irpnt + n
  if (ier < 0) return
  call split (accel,suba4,suba5,subq59,subq60,subq61,subq62,subq63,subq64, &
    subq65,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine mbic2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MBIC2 drives the block factorization (version 1, modified) method.
!
  external accel, suba1, subq25, copy, noadp
  external ibfcs3
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacs (2,ibfcs3,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  iwkpt1 = irpnt
  irpnt = irpnt + kblsz
  if (ier < 0) return
  call split (accel,suba1,suba1,subq25,subq25,subq25,subq25,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - kblsz
  return
end
subroutine mbic3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MBIC3 drives the block factorization (version 1, modified) method.
!
  external accel, suba4, suba5, subq70, subq71, subq72
  external subq73, subq74, subq75, noadp
  external ibfcn3
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacmz (2,ibfcn3,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*kblsz
  call split (accel,suba4,suba5,subq70,subq71,subq72,subq73,subq74,subq75, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*kblsz
  return
end
subroutine mbic7 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MBIC7 drives the block factorization (version 1, modified) method.
!
! (multi-color ordering)
!
  external accel, suba2, suba3, subq34, subq35, subq36
  external subq37, subq38, subq39, noadp
  external ibfcn3
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  t1 = timer (dummy)
  if (ifact == 1) call bfacmy (2,ibfcn3,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*ncmax
  call split (accel,suba2,suba3,subq34,subq35,subq36,subq37,subq38,subq39, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*ncmax
  return
end
subroutine mbicx2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MBICX2 drives the block factorization (version 2, modified) method.
!
  external accel, suba1, subq25, copy, noadp
  external ibfcs4
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacs (4,ibfcs4,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  iwkpt1 = irpnt
  irpnt = irpnt + kblsz
  if (ier < 0) return
  call split (accel,suba1,suba1,subq25,subq25,subq25,subq25,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - kblsz
  return
end
subroutine mbicx3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MBICX3 drives the block factorization (version 2, modified) method.
!
  external accel, suba4, suba5, subq70, subq71, subq72
  external subq73, subq74, subq75, noadp
  external ibfcn4
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  call blkdef (coef,jcoef,wksp,iwksp,n,ier)
  if (ier < 0) return
  t1 = timer (dummy)
  if (ifact == 1) call bfacmz (4,ibfcn4,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*kblsz
  call split (accel,suba4,suba5,subq70,subq71,subq72,subq73,subq74,subq75, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*kblsz
  return
end
subroutine mbicx7 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwkspiparm,rparm,ier)
!
!*******************************************************************************
!
!! MBICX7 drives the block factorization (version 2, modified method).
!
!  (multi-color ordering)
!
  external accel, suba2, suba3, subq34, subq35, subq36
  external subq37, subq38, subq39, noadp
  external ibfcn4
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  t1 = timer (dummy)
  if (ifact == 1) call bfacmy (4,ibfcn4,coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + 2*ncmax
  call split (accel,suba2,suba3,subq34,subq35,subq36,subq37,subq38,subq39, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - 2*ncmax
  return
end
subroutine mcopy (lda,ldb,n,m,a,b)
!
!*******************************************************************************
!
!! MCOPY copies an array into array.
!
!
!  Parameters:
!
!        lda     leading dimension of array a
!        ldb     leading dimension of array b
!        n       number of rows in a to be copied
!        m       number of columns in a to be copied
!        a,b     arrays
!
!  
!
  dimension a(lda,1), b(ldb,1)
!
  do 15 j = 1,m
     do 10 i = 1,n
 10      b(i,j) = a(i,j)
 15   continue
  return
end
subroutine me (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! ME is the user interface to the minimal error algorithm of Fridman.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call mew (suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs,wksp(irpnt), &
    nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine mew (suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs,wk,nw, &
  iparm,rparm,ier)
!
!*******************************************************************************
!
!! MEW runs the minimal error algorithm of Fridman.
!
!
! the reference is: v. m. fridman, "the method of minimum iterations
!.", ussr computational math. and math. phys., vol. 2, 1962,
! pp. 362-3.
!
! two-sided preconditioning is implemented.  the iteration matrix
! should be symmetric for this algorithm to work.
!
! we have introduced periodic scaling of the direction vectors, to
! prevent overflow.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subql, subqr
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! the following indexing functions are used to access the old
! direction vectors --
!
  indp(i)  = ip  + mod(i,2)*n
  indpt(i) = ipt + mod(i,2)*n
!
! various preliminary calculations.
!
  dot = 0e0
  nwusd = 0
  ier = 0
  iacel = 4
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  if (level >= 2) write (nout,496)
496   format (' me')
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
!
! initialize the stopping test.
!
  call inithv (0)
  zhave  = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr, coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx, wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
! memory allocation, etc.
!
! nomenclature -- r  -- residual of the original system.
!                 z  -- inv(ql)*r
!                 zt -- inv(qr)*z
!
  ip = 1
  ipt = ip + 2*n
  iz = ipt + 2*n
  ir = iz + n
  iv1 = ir + n
  if (.not. rcalp) iv1 = ir
  izt = iv1 + n
  iv2 = izt + n
  if (.not. ztcalp) iv2 = izt
  iqlap = iv1
  iqrlap = iv2
  iwfree = iv2 + n
!
! note that memory usage has been overlapped whenever possible,
! in order to save space.
!
  nwusd = max(nwusd,iwfree-1)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  rhave = rcalp
  zthave = ztcalp
!
! perform first-iterate calculations
!
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(ir))
  call vexopy (n,wk(ir),rhs,wk(ir),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(ir),wk(iz))
  call subqr (coef,jcoef,wfac,jwfac,n,wk(iz),wk(izt))
!
!  Begin iteration loop.
!
! determine whether or not to stop --
! note that we have already done the calculations necessary so that suba
! and subql are not actually used by pstop.
!
 10   call inithv (1)
  nwpstp = nw - (iwfree-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    wk(ir),wk(iz),wk(izt),wk(iwfree),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iwfree-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
! compute p(n), the direction vector, and inv(qr)*p(n) (=pt(n)).
!
  scal = 1e0
!
! first, case of in == 0
!
  if (in /= 0) go to 100
  toplam = vdot (n,wk(iz),wk(iz))
  call suba (coef,jcoef,wfac,jwfac,n,wk(izt),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(indp(in)))
  call subqr (coef,jcoef,wfac,jwfac,n,wk(indp(in)),wk(indpt(in)))
  go to 120
!
! case in > 0
!
 100  toplam = vdot (n,wk(indp(in-1)),wk(iz))
  bet1 = - vdot (n,wk(indp(in-1)),wk(iqlap)) / dot
  if (in /= 1) go to 110
!
! case in == 1
!
  call vtriad (n,wk(indp(in)),wk(iqlap),bet1,wk(indp(in-1)),1)
  call vtriad (n,wk(indpt(in)),wk(iqrlap),bet1,wk(indpt(in-1)),1)
  go to 120
!
! case in > 1
!
 110  bet2 = - vdot (n,wk(indp(in-2)),wk(iqlap)) / dotold
  call vtriad (n,wk(indp(in)), wk(iqlap), bet2,wk(indp(in-2)), 1)
  call vtriad (n,wk(indpt(in)),wk(iqrlap),bet2,wk(indpt(in-2)),1)
  call vtriad (n,wk(indp(in)), wk(indp(in)), bet1,wk(indp(in-1)), 1)
  call vtriad (n,wk(indpt(in)),wk(indpt(in)),bet1,wk(indpt(in-1)),1)
!
! at this point, we are finished forming the latest direction vector.
! we proceed to calculate lambda and update the solution and the
! residual.
!
 120  dotold = dot
  dot = vdot (n,wk(indp(in)),wk(indp(in)))
!     if (dot < srelpr) go to 998
!
! scale direction vector if necessary.
  if (dot<srelpr**2 .or. dot>(1e0/srelpr**2)) then
    scal = sqrt(dot)
    call vtriad (n,wk(indp(in)), xxx,1e0/scal,wk(indp(in)), 2)
    call vtriad (n,wk(indpt(in)),xxx,1e0/scal,wk(indpt(in)),2)
    dot = 1e0
  end if
!
 124  vlamda = toplam / dot / scal
!
! u --
!
  call vtriad (n,u,u,vlamda,wk(indpt(in)),1)
!
! r --
!
  call suba (coef,jcoef,wfac,jwfac,n,wk(indpt(in)),wk(iv2))
  if (rhave) call vtriad (n,wk(ir),wk(ir),-vlamda,wk(iv2),1)
!
! z --
!
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv2),wk(iqlap))
  call vtriad (n,wk(iz),wk(iz),-vlamda,wk(iqlap),1)
!
! zt --
!
  call subqr (coef,jcoef,wfac,jwfac,n,wk(iqlap),wk(iqrlap))
  if (zthave) call vtriad (n,wk(izt),wk(izt),-vlamda,wk(iqrlap),1)
!
! proceed to next iteration
!
  in = in + 1
  is = is + 1
  go to 10
!
!  Finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'mew')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' me converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 997  call ershow (ier,'mew')
  go to 735
!
 998  ier = -15
  call ershow (ier,'mew')
  go to 725
!
 999  ier = -2
  call ershow (ier,'mew')
  go to 735
!
end
subroutine mfact (coef,jcoef,wksp,iwksp,nn,ier)
!
!*******************************************************************************
!
!! MFACT computes a line factorization of a multi-color matrix.
!
!
!  Parameters:
!
!        n        problem size
!        nfactr   factorization size
!
!  common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
!
!  check for sufficient workspace to store factor.
!
  n = nn
  ndt = iwksp(iblock+2) - 1
  ndb = iwksp(iblock+ncolor*3+2)
  nwdiag = ndt + ndb + 1
  nfactr = n*nwdiag
  call needw ('mfact',0,irpnt,nfactr,ier)
  if (ier < 0) return
!
  ifactr = irpnt
  do 15 j = 1,nwdiag
     ipt1 = (j - 1)*ndim + 1
     ipt2 = (j - 1)*n + ifactr
     call vcopy (n,coef(ipt1),wksp(ipt2))
 15   continue
!
!  factor.
!
  call bdfac (n,n,n,ndt,ndb,wksp(ifactr),1)
  irpnt = irpnt + nfactr
  return
end
subroutine mic1 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MIC1 drives the MIC preconditioner.
!
  external accel, suba8, suba9, subq86, subq87, subq88
  external subq89, subq90, subq91, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
  n = nn
  if (ifact == 0 .and. lvfill > 0) go to 20
  call move1 (ndim,mdim,n,maxnz,jcoef,coef,maxt,maxb,ier)
  if (ier < 0) then
     call ershow (ier,'mic1')
     return
  end if
 20   t1 = timer (dummy)
  if (ifact == 1) call pfact1 (coef,jcoef,wksp,iwksp,n,2,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba8,suba9,subq86,subq87,subq88,subq89,subq90,subq91, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine mic2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MIC2 drives the symmetric MIC preconditioner.
!
  external accel, suba1, subq13, subq14, subq15, subq16, subq17, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
!
  t1 = timer (dummy)
  if (ifact == 1) call pfact2 (coef,jcoef,wksp,iwksp,n,2,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  leniw = max (maxnz,nfacti)
  iwkpt1 = iipnt
  iipnt = iipnt + leniw
  call split (accel,suba1,suba1,subq13,subq13,subq14,subq15,subq16,subq17, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  iipnt = iipnt - leniw
  return
end
subroutine mic3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MIC3 drives the nonsymmetric MIC preconditioner.
!
  external accel, suba4, suba5, subq48, subq49, subq50
  external subq51, subq52, subq53, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
!
  n = nn
  call needw ('mic3',1,iipnt,maxnz,ier)
  if (ier < 0) return
  call needw ('mic3',0,irpnt,n,ier)
  if (ier < 0) return
  if (ifact == 0 .and. lvfill > 0) go to 20
  call move2 (ndim,n,maxnz,jcoef,coef,wksp(irpnt),iwksp(iipnt),maxt,maxb)
 20   t1 = timer (dummy)
  if (ifact == 1) call pfact3 (coef,jcoef,wksp,iwksp,n,2,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  leniw = max (maxnz,nfacti)
  iwkpt1 = iipnt
  iipnt = iipnt + leniw
  call split (accel,suba4,suba5,subq48,subq49,subq50,subq51,subq52,subq53, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  iipnt = iipnt - leniw
  return
end
subroutine mic6 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! MIC6 drives the MIC preconditioner.
!
!     (multi-color ordering)
!
  external accel, suba8, suba9, sub104, sub105, sub106
  external sub107, sub108, sub109, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
  n = nn
  t1 = timer (dummy)
  if (ifact == 1) call pfactc (coef,jcoef,wksp,iwksp,n,2,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba8,suba9,sub104,sub105,sub106,sub107,sub108,sub109, &
    noadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine move1 (ndim,mdim,nn,maxnzz,jcoef,coef,nt,nb,ier)
!
!*******************************************************************************
!
!! MOVE1 moves the data structure to the form d/t/b.
!
!
!     d is the main diagonal, the t columns contain only upper
!     triangular elements and the b columns contain only lower
!     triangular elements.  thus the upper and lower triangle
!     elements are segregated into separate columns of coef,
!     with the upper elements coming first.
!     (Purdue data structure, natural ordering, with point
!     ic or point ssor preconditionings)
!
!  Parameters:
!
!         ndim     row dimension of coef array in defining routine
!         mdim     column dimension of coef array in defining routine
!         n        order of system (= nn)
!         maxnz    number of columns in coef array (= maxnzz)
!         jcoef    integer matrix representation array
!         coef     matrix representation array
!         nt       number of columns needed to store t, the upper
!                   triangular part of coef
!         nb       number of columns needed to store b, the lower
!                   triangular part of coef
!         ier      error code
!                  =  0  no errors detected
!                  = -9  mdim < 1+nt+nb.  hence insufficient room
!                        to store adjusted matrix
!
!  
!
  integer   jcoef(ndim,1)
  dimension coef(ndim,1)

  n = nn
  maxnz = maxnzz
!
!  determine maximum number of nonzeros per row in t and b.
!
  ntt = 0
  nbb = 0
  if (maxnz <= 1) go to 999
  do 25 i = 1,n
     ntrow = 0
     nbrow = 0
     do 20 j = 2,maxnz
        if (jcoef(i,j) - i) 10,20,15
 10         nbrow = nbrow + 1
        go to 20
 15         ntrow = ntrow + 1
 20      continue
     if (ntrow > ntt) ntt = ntrow
     if (nbrow > nbb) nbb = nbrow
 25   continue
!
!  shuffle matrix so that t is first.
!
  ndtb = ntt + nbb + 1
  if (ndtb <= mdim) go to 30
!
!  error -- mdim is too small.
!
  ier = -9
  go to 999
!
!  permute elements of each row.
!
 30   if (ntt*nbb == 0) go to 999
  if (ndtb <= maxnz) go to 40
  maxz = maxnz + 1
  do 35 j = maxz,ndtb
  do 35 i = 1,n
     coef(i,j) = 0.0
     jcoef(i,j) = i
 35   continue
  maxnz = ndtb
 40   nt2 = ntt + 1
  nb1 = nt2 + 1
  do 65 i = 1,n
     jbc = nt2
     do 50 jtc = 2,nt2
        if (jcoef(i,jtc) >= i) go to 50
 45         jbc = jbc + 1
        if (jcoef(i,jbc) < i) go to 45
        jtemp = jcoef(i,jtc)
        jcoef(i,jtc) = jcoef(i,jbc)
        jcoef(i,jbc) = jtemp
        temp = coef(i,jtc)
        coef(i,jtc) = coef(i,jbc)
        coef(i,jbc) = temp
 50      continue
     jtc = 1
     do 60 jbc = nb1,maxnz
        if (jcoef(i,jbc) <= i) go to 60
 55         jtc = jtc + 1
        if (jcoef(i,jtc) > i) go to 55
        jtemp = jcoef(i,jtc)
        jcoef(i,jtc) = jcoef(i,jbc)
        jcoef(i,jbc) = jtemp
        temp = coef(i,jtc)
        coef(i,jtc) = coef(i,jbc)
        coef(i,jbc) = temp
 60      continue
 65   continue
!
!  exit.
!
 999  nt = ntt
  nb = nbb
  maxnzz = maxnz
  return
end
subroutine move2 (ndim,nn,maxnzz,jcoef,coef,work,iwork,nt,nb)
!
!*******************************************************************************
!
!! MOVE2 moves the data structure to the form d/t/b.
!
!
!     d is the main diagonal, the t columns contain only upper
!     triangular elements and the b columns contain only lower
!     triangular elements.  thus the upper and lower triangle
!     elements are segregated into separate columns of coef,
!     with the upper elements coming first.
!     (diagonal data structure, natural ordering, with point
!     ic or point ssor preconditionings)
!
!  Parameters:
!
!         ndim     row dimension of coef array in defining routine
!         n        order of system (= nn)
!         maxnz    number of columns in coef array (= maxnzz)
!         jcoef    integer matrix representation array
!         coef     matrix representation array
!         work     real workspace array of length n
!         iwork    integer work array of length maxnz
!         nt       number of columns needed to store t, the upper
!                   triangular part of coef
!         nb       number of columns needed to store b, the lower
!                   triangular part of coef
!
!  
!
  integer   jcoef(2), iwork(1)
  dimension coef(ndim,1), work(1)
!
  n = nn
  maxnz = maxnzz
  ntt = 0
  nbb = 0
  if (maxnz <= 1) go to 999
!
!  compute nbb and ntt.
!
  do 10 j = 1,maxnz
     ndiag = jcoef(j)
     if (ndiag > 0) ntt = ntt + 1
     if (ndiag < 0) nbb = nbb + 1
 10   continue
!
!  compute pointers into sorted jcoef.
!
!  code jcoef.
!
  do 15 j = 1,maxnz
     if (jcoef(j) < 0) jcoef(j) = n - jcoef(j)
 15   continue
  iwork(1) = 1
  do 30 j = 2,maxnz
     iaux = jcoef(j)
     do 20 k = 1,j-1
        i = j - k
        ktemp = iwork(i)
        if (iaux > jcoef(ktemp)) go to 25
        iwork(i+1) = iwork(i)
 20      continue
     i = 0
 25      iwork(i+1) = j
 30   continue
!
!  decode jcoef.
!
  do 35 j = 1,maxnz
     if (jcoef(j) > n) jcoef(j) = n - jcoef(j)
 35   continue
!
!  sort coef and jcoef.
!
  do 40 i = 1,maxnz
     if (iwork(i) == i) iwork(i) = 0
 40   continue
  do 65 ii = 1,maxnz
     k = iwork(ii)
     if (k == 0) go to 65
     i = ii
 45      jtemp = jcoef(i)
     jcoef(i) = jcoef(k)
     jcoef(k) = jtemp
     do 50 l = 1,n
        work(l) = coef(l,i)
        coef(l,i) = coef(l,k)
        coef(l,k) = work(l)
 50      continue
     iwork(i) = 0
     do 55 j = ii,maxnz
        if (iwork(j) == i) go to 60
 55      continue
     go to 65
 60      i = j
     if (i /= k) go to 45
     iwork(k) = 0
 65   continue
!
!  exit.
!
 999  nt = ntt
  nb = nbb
  return
end
subroutine move3 (ndim,mdim,nn,maxnzz,jcoef,coef,nt,nb,ncolor,nc,ier)
!
!*******************************************************************************
!
!! MOVE3 moves the data structure to the form d/t/b.
!
!
!     d is the main diagonal, the t columns contain only upper
!     triangular elements and the b columns contain only lower
!     triangular elements.  thus the upper and lower triangle
!     elements are segregated into separate columns of coef,
!     with the upper elements coming first.
!     the above segregation is done for each color.
!     (Purdue data structure, multi-color ordering, with point
!     ic or point ssor preconditionings)
!
!  Parameters:
!
!         ndim     row dimension of coef array in defining routine
!         mdim     column dimension of coef array in defining routine
!         n        order of system (= nn)
!         maxnz    number of columns in coef array (= maxnzz)
!         jcoef    integer matrix representation array
!         coef     matrix representation array
!         nt       integer vector of length ncolor.  for each color,
!                   the number of columns needed to store t, the upper
!                   triangular part of the matrix for those rows.
!         nb       integer vector of length ncolor.  for each color,
!                   the number of columns needed to store b, the lower
!                   triangular part of the matrix for those rows.
!         ncolor   number of colors
!         nc       integer vector of length ncolor, giving the number
!                   of nodes for each color.
!         ier      error code
!                  =  0  no errors detected
!                  = -9  mdim < 1+nt+nb.  hence insufficient room
!                        to store adjusted matrix
!
!  
!
  integer   jcoef(ndim,1), nt(1), nb(1), nc(1)
  dimension coef(ndim,1)
!
!
  n = nn
  maxnz = maxnzz
!
  ist = 1
  do 85 icol = 1,ncolor
     ncol = nc(icol)
     ied = ist + ncol - 1
!
!  determine maximum number of nonzeros per row in t and b.
!
     ntt = 0
     nbb = 0
     if (maxnz <= 1) go to 80
     do 25 i = ist,ied
        ntrow = 0
        nbrow = 0
        do 20 j = 2,maxnz
           if (jcoef(i,j) - i) 10,20,15
 10            nbrow = nbrow + 1
           go to 20
 15            ntrow = ntrow + 1
 20         continue
        if (ntrow > ntt) ntt = ntrow
        if (nbrow > nbb) nbb = nbrow
 25      continue
!
!  shuffle matrix so that t is first.
!
     ndtb = ntt + nbb + 1
     if (ndtb <= mdim) go to 30
!
!  error -- mdim is too small.
!
     ier = -9
     go to 999
!
!  permute elements of each row.
!
 30      if (ndtb <= maxnz) go to 40
     maxz = maxnz + 1
     do 35 j = maxz,ndtb
     do 35 i = 1,n
        coef(i,j) = 0.0
        jcoef(i,j) = i
 35      continue
     maxnz = ndtb
 40      nt2 = ntt + 1
     nb1 = ntt + 2
     nz1 = 2 + ntt + nbb
     do 75 i = ist,ied
        jbc = nt2
        do 50 jtc = 2,nt2
           if (jtc > nt2) go to 50
           if (jcoef(i,jtc) >= i) go to 50
 45            jbc = jbc + 1
           if (jcoef(i,jbc) < i) go to 45
           jtemp = jcoef(i,jtc)
           jcoef(i,jtc) = jcoef(i,jbc)
           jcoef(i,jbc) = jtemp
           temp = coef(i,jtc)
           coef(i,jtc) = coef(i,jbc)
           coef(i,jbc) = temp
 50         continue
        jtc = 1
        do 60 jbc = nb1,maxnz
           if (jbc > maxnz) go to 60
           if (jcoef(i,jbc) <= i) go to 60
 55            jtc = jtc + 1
           if (jcoef(i,jtc) > i) go to 55
           jtemp = jcoef(i,jtc)
           jcoef(i,jtc) = jcoef(i,jbc)
           jcoef(i,jbc) = jtemp
           temp = coef(i,jtc)
           coef(i,jtc) = coef(i,jbc)
           coef(i,jbc) = temp
 60         continue
        jbc = nt2
        do 70 jzc = nz1,maxnz
           if (jzc > maxnz) go to 70
           if (jcoef(i,jzc) >= i) go to 70
 65            jbc = jbc + 1
           if (jcoef(i,jbc) < i) go to 65
           jtemp = jcoef(i,jzc)
           jcoef(i,jzc) = jcoef(i,jbc)
           jcoef(i,jbc) = jtemp
           temp = coef(i,jzc)
           coef(i,jzc) = coef(i,jbc)
           coef(i,jbc) = temp
 70         continue
 75      continue
!
 80      nt(icol) = ntt
     nb(icol) = nbb
     ist = ist + ncol
 85   continue
!
!  exit.
!
 999  maxnzz = maxnz
  return
end
subroutine move4 (ndim,nn,maxnew,jcnew,coef,ncol,nc,work,iwork)
!
!*******************************************************************************
!
!! MOVE4 moves the data structure to the form dc/tc/bc.
!
!
!     dc is the main diagonal block, tc is the upper triangular
!     block matrices, and db is the lower triangular block
!     matrices.
!     the above segregation is done for each color.
!     (diagonal data structure, multi-color ordering, with
!     ic or ssor preconditionings (point or block))
!
!  Parameters:
!
!         ndim     row dimension of coef array in defining routine
!         n        order of system (= nn)
!         maxnew   integer vector giving the number of diagonals
!                   created for each color
!         jcnew    integer array of size ncolor*max(maxnew(i))
!                   giving the diagonal numbers for each color
!         coef     matrix representation array
!         ncolor   number of colors
!         nc       integer vector of length ncolor, giving the number
!                   of nodes for each color.
!         work     real workspace array of length max (nc(i))
!         iwork    integer work array of length max (maxnew(i))
!
!  
!
  integer   maxnew(1), jcnew(ncol,1), nc(1), iwork(1)
  dimension coef(ndim,1), work(1)
!
  n = nn
  ncolor = ncol
  ist = 1
  do 70 icol = 1,ncolor
     ncc = nc(icol)
     ied = ist + ncc - 1
!
!  compute pointers into sorted jcnew.
!
!  code jcnew.
!
     maxnz = maxnew(icol)
     do 15 j = 1,maxnz
        do 5 i = ist,ied
           if (coef(i,j) /= 0.0) go to 10
 5          continue
        go to 15
 10         jd = jcnew(icol,j)
        jcol = i + jd
        if (jcol < i .and. jcol >= ist) jcnew(icol,j) = n - jd
        if (jcol > ied) jcnew(icol,j) = 2*n + jd
        if (jcol < ist) jcnew(icol,j) = 3*n - jd
 15      continue
     iwork(1) = 1
     do 30 j = 2,maxnz
        iaux = jcnew(icol,j)
        do 20 k = 1,j-1
           i = j - k
           ktemp = iwork(i)
           if (iaux > jcnew(icol,ktemp)) go to 25
           iwork(i+1) = iwork(i)
 20         continue
        i = 0
 25         iwork(i+1) = j
 30      continue
!
!  decode jcnew.
!
     do 35 j = 1,maxnz
        jd = jcnew(icol,j)
        if (jd > n .and. jd < 2*n) jcnew(icol,j) = n - jd
        if (jd > 2*n .and. jd < 3*n) jcnew(icol,j) = jd - 2*n
        if (jd > 3*n) jcnew(icol,j) = 3*n - jd
 35      continue
!
!  sort coef and jcnew.
!
     do 40 i = 1,maxnz
        if (iwork(i) == i) iwork(i) = 0
 40      continue
     do 65 ii = 1,maxnz
        k = iwork(ii)
        if (k == 0) go to 65
        i = ii
 45         jtemp = jcnew(icol,i)
        jcnew(icol,i) = jcnew(icol,k)
        jcnew(icol,k) = jtemp
        do 50 l = ist,ied
           work(l-ist+1) = coef(l,i)
           coef(l,i) = coef(l,k)
           coef(l,k) = work(l-ist+1)
 50         continue
        iwork(i) = 0
        do 55 j = ii,maxnz
           if (iwork(j) == i) go to 60
 55         continue
        go to 65
 60         i = j
        if (i /= k) go to 45
        iwork(k) = 0
 65      continue
     ist = ist + ncc
 70   continue
!
!  exit.
!
  return
end
subroutine move5 (ndim,n,maxnz,jcoef,coef)
!
!*******************************************************************************
!
!! MOVE5 moves the data structure to the form dc/tc/bc.
!
!
!     dc is the main diagonal block, tc is the upper triangular
!     block matrices, and db is the lower triangular block
!     matrices.
!     (diagonal data structure, with constant block size)
!
!  Parameters:
!
!         ndim     row dimension of coef array in defining routine
!         n        order of system
!         maxnz    number of diagonals stored
!         jcoef    integer vector of length maxnz giving the
!                   diagonal numbers
!         coef     matrix representation array
!
!  
!
  dimension coef(ndim,maxnz), jcoef(maxnz)
!
!  move dc to the first columns.
!
  jsh = 1
  jcol = 1
  jget = 0
 5    do 10 j = 1,maxnz
     jd = jcoef(j)
     if (jd == jget) go to 15
 10   continue
  if (jsh < 0) go to 30
  jsh = -1
  jget = -1
  go to 5
 15   if (j == jcol) go to 25
  do 20 i = 1,n
     temp = coef(i,j)
     coef(i,j) = coef(i,jcol)
     coef(i,jcol) = temp
 20   continue
  jcoef(j) = jcoef(jcol)
  jcoef(jcol) = jd
 25   jcol = jcol + 1
  jget = jget + jsh
  go to 5
!
!  move tc, bc to the next columns.
!
 30   if (jcol > maxnz) return
  do 35 j = jcol,maxnz
     jd = jcoef(j)
     if (jd < 0) jcoef(j) = n - jd
 35   continue
  jcolsv = jcol
 40   jsml = jcol
  do 45 j = jcol,maxnz
     jd = jcoef(j)
     if (jd < jcoef(jsml)) jsml = j
 45   continue
  if (jsml == jcol) go to 55
  do 50 i = 1,n
     temp = coef(i,jsml)
     coef(i,jsml) = coef(i,jcol)
     coef(i,jcol) = temp
 50   continue
  jtemp = jcoef(jsml)
  jcoef(jsml) = jcoef(jcol)
  jcoef(jcol) = jtemp
 55   jcol = jcol + 1
  if (jcol <= maxnz) go to 40
!
!  uncode jcoef.
!
  do j = jcolsv,maxnz
    jd = jcoef(j)
    if (jd > n) jcoef(j) = n - jd
  end do

  return
end
subroutine mul1t (ndim,maxnz,coef,jcoef,wksp,nn,x,y)
!
!*******************************************************************************
!
!! MUL1T computes y = (A**t)*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in column one.
!     (Purdue storage format)
!
!  Parameters:
!
!         ndim     row dimension of coef in defining routine
!         maxnz    number of columns in coef
!         coef     array of matrix nonzeros
!         jcoef    array of matrix column numbers
!         wksp     workspace array of length n
!         n        dimension of matrix (= nn)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension coef(ndim,2), x(1), y(1), wksp(1)
  integer   jcoef(ndim,2)
!
  n = nn

  do i = 1,n
    y(i) = coef(i,1)*x(i)
  end do

  if (maxnz <= 1) return
  maxm1 = maxnz - 1
  call vaddpt (ndim,ndim,n,maxm1,coef(1,2),jcoef(1,2),y,x,wksp)
  return
end
subroutine mul2nt (ndim,maxnz,coef,jcoef,nn,x,y)
!
!*******************************************************************************
!
!! MUL2NT computes y = (A**t)*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in column one.  all diagonals of
!     the matrix must be stored.
!     (nonsymmetric diagonal storage format)
!
!  Parameters:
!
!         ndim     row dimension of coef in defining routine
!         maxnz    number of columns in coef
!         coef     array of matrix diagonals
!         jcoef    array of matrix diagonal numbers
!         n        dimension of matrix (= nn)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension coef(ndim,2), x(1), y(1)
  integer   jcoef(2)
!
  n = nn

  do i = 1,n
    y(i) = coef(i,1)*x(i)
  end do

  if (maxnz <= 1) return
  maxm1 = maxnz - 1
  call vadddt (ndim,1,n,n,maxm1,coef(1,2),jcoef(2),y,x,0)
  return
end
subroutine mul3nt (mm,np,a,ia,ja,wksp,x,y)
!
!*******************************************************************************
!
!! MUL3NT computes y = (A**t)*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in the first partition.
!     (nonsymmetric sparse storage format)
!
!  Parameters:
!
!         m        number of partitions
!         np       integer vector of length m+1 giving partition
!                    pointers
!         a        real vector giving matrix coefficients
!         ia       integer vector giving i values
!         ja       integer vector giving j values
!         wksp     workspace vector of length 2*n (keygs = 1 only)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension a(1), x(1), y(1), wksp(1)
  integer   np(2), ia(1), ja(1)
!
  m = mm
  ied = np(2) - 1

  do i = 1,ied
    y(i) = a(i)*x(i)
  end do

  mm1 = m - 1
  call vadds (mm1,np(2),ja,ia,a,y,x,wksp)
  return
end
subroutine muldc (ndim,nn,coef,ncolor,nc,maxnew,jcnew,x,y)
!
!*******************************************************************************
!
!! MULDC computes y = A*x for a matrix permuted to an ncolor x ncolor block matrix..
!
!
!  The matrix is stored in diagonal format.
!
!  Parameters:
!
!        ndim      row dimension of coef array
!        n         order of system
!        coef      real array of coefficients
!        ncolor    number of colors in the permutation (= ncol)
!        nc        integer vector of length ncolor giving the
!                   number of nodes for each color
!        maxnew    integer vector giving the number of diagonals
!                   created for each color
!        jcnew     integer array of size ncolor*max(maxnew(i))
!                   giving the diagonal numbers for each color
!        x         vector of length n to be multiplied by
!        y         vector of length n to contain result vector
!
!  
!
  integer nc(1), maxnew(1), jcnew(ncolor,2)
  dimension coef(ndim,2), x(1), y(1)
!
  n = nn

  do i =1,n
    y(i) = coef(i,1)*x(i)
  end do

  i1 = 1
  joff = 0
  do 15 k = 1,ncolor
     ncc = nc(k)
     jlim = maxnew(k) - 1
     call vaddd (ndim,ncolor,ncc,n,jlim,coef(i1,2),jcnew(k,2), y(i1),x,joff)
     i1 = i1 + ncc
     joff = joff - ncc
 15   continue
  return
end
subroutine muldct (ndim,nn,coef,ncolor,nc,maxnew,jcnew,x,y)
!
!*******************************************************************************
!
!! MULDCT computes y = (A**t)*x for a matrix permuted to an ncolor x ncolor block matrix.
!
!
!  The matrix is stored in diagonal format.
!
!  Parameters:
!
!        ndim      row dimension of coef array
!        n         order of system
!        coef      real array of coefficients
!        ncolor    number of colors in the permutation (= ncol)
!        nc        integer vector of length ncolor giving the
!                   number of nodes for each color
!        maxnew    integer vector giving the number of diagonals
!                   created for each color
!        jcnew     integer array of size ncolor*max(maxnew(i))
!                   giving the diagonal numbers for each color
!        x         vector of length n to be multiplied by
!        y         vector of length n to contain result vector
!
!  
!
  integer nc(1), maxnew(1), jcnew(ncolor,2)
  dimension coef(ndim,2), x(1), y(1)
!
  n = nn
  do i =1,n
    y(i) = coef(i,1)*x(i)
  end do

  i1 = 1
  joff = 0
  do 15 k = 1,ncolor
     ncc = nc(k)
     jlim = maxnew(k) - 1
     call vadddt (ndim,ncolor,ncc,n,jlim,coef(i1,2),jcnew(k,2),y,x(i1),joff)
     i1 = i1 + ncc
     joff = joff - ncc
 15   continue
  return
end
subroutine mult1 (ndim,maxnz,coef,jcoef,wksp,nn,x,y)
!
!*******************************************************************************
!
!! MULT1 computes y = A*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in column one.
!     (Purdue storage format)
!
!  Parameters:
!
!         ndim     row dimension of coef in defining routine
!         maxnz    number of columns in coef
!         coef     array of matrix nonzeros
!         jcoef    array of matrix column numbers
!         wksp     workspace array of length n
!         n        order of matrix (= nn)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension coef(ndim,2), x(1), y(1), wksp(1)
  integer   jcoef(ndim,2)
!
  n = nn
  maxm1 = maxnz - 1
  do 10 i = 1,n
 10   y(i) = coef(i,1)*x(i)
  call vaddp (ndim,ndim,n,maxm1,coef(1,2),jcoef(1,2),y,x,wksp)
  return
end
subroutine mult2n (ndim,maxnz,coef,jcoef,nn,x,y)
!
!*******************************************************************************
!
!! MULT2N computes y = A*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in column one.  all diagonals of
!     the matrix must be stored.
!     (nonsymmetric diagonal storage format)
!
!  Parameters:
!
!         ndim     row dimension of coef in defining routine
!         maxnz    number of columns in coef
!         coef     array of matrix diagonals
!         jcoef    array of matrix diagonal numbers
!         n        dimension of matrix (= nn)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension coef(ndim,2), x(1), y(1)
  integer   jcoef(2)
!
  n = nn
  do i = 1,n
    y(i) = coef(i,1)*x(i)
  end do

  if (maxnz <= 1) return
  maxm1 = maxnz - 1
  call vaddd (ndim,1,n,n,maxm1,coef(1,2),jcoef(2),y,x,0)
  return
end
subroutine mult2s (ndim,maxnz,coef,jcoef,nn,x,y)
!
!*******************************************************************************
!
!! MULT2S computes y = A*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in column 1.  only the upper
!     diagonals and the main diagonal are assumed stored.
!     (symmetric diagonal storage format)
!
!  Parameters:
!
!         ndim     row dimension of coef in defining routine
!         maxnz    number of columns in coef
!         coef     array of matrix diagonals
!         jcoef    array of matrix diagonal numbers
!         n        dimension of matrix (= nn)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension coef(ndim,1), x(1), y(1)
  integer   jcoef(2)
!
  n = nn

  do i = 1,n
    y(i) = coef(i,1)*x(i)
  end do

  if (maxnz <= 1) return

  do 25 j = 2,maxnz
     ind = jcoef(j)
     len = n - ind
     do 15 i = 1,len
 15      y(i) = y(i) + coef(i,j)*x(i+ind)
     do 20 i = 1,len
 20      y(i+ind) = y(i+ind) + coef(i,j)*x(i)
 25   continue
  return
end
subroutine mult3 (mm,np,a,ia,ja,wksp,x,y)
!
!*******************************************************************************
!
!! MULT3 computes y = A*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in the first partition.
!     (symmetric sparse storage format)
!
!  Parameters:
!
!         m        number of partitions
!         np       integer vector of length m+1 giving partition
!                    pointers
!         a        real vector giving matrix coefficients
!         ia       integer vector giving i values
!         ja       integer vector giving j values
!         wksp     workspace vector of length 2*n (keygs = 1 only)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension a(1), x(1), y(1), wksp(1)
  integer   np(2), ia(1), ja(1)
!
  m = mm
  ied = np(2) - 1
  do i = 1, ied
    y(i) = a(i)*x(i)
  end do

  mm1 = m - 1
  call vadds (mm1,np(2),ia,ja,a,y,x,wksp)
  call vadds (mm1,np(2),ja,ia,a,y,x,wksp)

  return
end
subroutine mult3n (mm,np,a,ia,ja,wksp,x,y)
!
!*******************************************************************************
!
!! MULT3N computes y = A*x, a matrix-vector product.
!
!
!     the diagonal is assumed to be in the first partition.
!     (nonsymmetric sparse storage format)
!
!  Parameters:
!
!         m        number of partitions
!         np       integer vector of length m+1 giving partition
!                    pointers
!         a        real vector giving matrix coefficients
!         ia       integer vector giving i values
!         ja       integer vector giving j values
!         wksp     workspace vector of length 2*n (keygs = 1 only)
!         x        multiplying vector of length n
!         y        product vector of length n
!
!  
!
  dimension a(1), x(1), y(1), wksp(1)
  integer   np(2), ia(1), ja(1)
!
  m = mm
  ied = np(2) - 1

  do i = 1,ied
    y(i) = a(i) * x(i)
  end do

  mm1 = m - 1
  call vadds (mm1,np(2),ia,ja,a,y,x,wksp)

  return
end
subroutine needw ( subnam, isw, istart, length, ier )
!
!*******************************************************************************
!
!! NEEDW determines if enough integer or real workspace is available.
!
!
!  Parameters:
!
!        subnam  name of calling routine
!        isw     switch for real or integer workspace check
!                 = 0     real
!                 = 1     integer
!        istart  starting address
!        length  length desired
!        ier     error indicator (output)
!                 = -2    insufficient real workspace
!                 = -3    insufficient integer workspace
!
!  
!
  character*10 subnam
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  newlen = istart + length - 1

  if ( isw /= 1 ) then

    if ( lenr < newlen ) then
      write ( *, * ) ' '
      write ( *, * ) 'NEEDW - Insufficient real workspace.'
      write ( *, * ) '  Increase needed is ', newlen - lenr
      ier = -2
      call ershow ( ier, subnam )
    end if

    irmax = max ( irmax, newlen )

  else

    if ( leni < newlen ) then
      write ( *, * ) ' '
      write ( *, * ) 'NEEDW - Insufficient integer workspace.'
      write ( *, * ) '  Increase needed is ', newlen - leni
      ier = -3
      call ershow ( ier, subnam )
    end if

    iimax = max ( iimax, newlen )

  end if

  return
end
subroutine neu1 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! NEU1 drives the Neumann polynomial preconditioner.
!
  external accel, suba8, suba9, subq94, subq95, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom4 / keygs, srelpr, keyzer
!
  n = nn
  call needw ('neu1',0,irpnt,n,ier)
  if (ier < 0) return
  iwkpt2 = irpnt
  irpnt = irpnt + n
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + n
  call split (accel,suba8,suba9,subq94,subq95,subq94,subq95,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  if (keygs == 1) irpnt = irpnt - n
  return
end
subroutine neu2 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! NEU2 drives the Neumann polynomial preconditioner.
!
  external accel, suba1, subq19, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  call needw ('neu2',0,irpnt,n,ier)
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba1,suba1,subq19,subq19,subq19,subq19,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine neu3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! NEU3 drives the Neumann polynomial preconditioner.
!
  external accel, suba4, suba5, subq56, subq57, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  call needw ('neu3',0,irpnt,n,ier)
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba4,suba5,subq56,subq57,subq56,subq57,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine neu4 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! NEU4 drives the Neumann polynomial preconditioner.
!
  external accel, suba12, sub111, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom4 / keygs, srelpr, keyzer
!
  n = nn
  call needw ('neu4',0,irpnt,n,ier)
  if (ier < 0) return
  iwkpt2 = irpnt
  irpnt = irpnt + n
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba12,suba12,sub111,sub111,sub111,sub111,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine neu5 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! NEU5 drives the Neumann polynomial preconditioner.
!
  external accel, suba13, suba14, sub114, sub115, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom4 / keygs, srelpr, keyzer
!
  n = nn
  call needw ('neu5',0,irpnt,n,ier)
  if (ier < 0) return
  iwkpt2 = irpnt
  irpnt = irpnt + n
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba13,suba14,sub114,sub115,sub114,sub115,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine nmcalc (coef,jcoef,wfac,jwfac,icall,subq,nn,rhs,ubar,wksp,ier)
!
!*******************************************************************************
!
!! NMCALC calculates the quantities
!
!       bnorm   = sqrt (rhs,rhs)
!       bnorm1  = any other norm of rhs needed for the stopping test
!       ubarnm  = sqrt (ubar,ubar)
!
!     which are needed in the stopping tests.
!
!     the stopping tests are --
!
!    (1)  (emax/emin) * sqrt ( (r ,zt)/(rhs,inv(q)*rhs) )
!    (2)  ( 1.0/emin) * sqrt ( (zt,zt)/(u,u) )
!    (3)  (emax/emin) * sqrt ( (zt,zt)/(inv(q)*rhs,inv(q)*rhs) )
!    (4)                sqrt ( (zt,zt)/(inv(q)*rhs,inv(q)*rhs) )
!    (5)                sqrt ( (r ,r )/(rhs,rhs) )
!    (6)                sqrt ( (u-ubar,u-ubar)/(ubar,ubar) )
!    (7)  (emax/emin) * sqrt ( (r,z)/(rhs,inv(ql)*rhs) )
!    (8)  ( 1.0/emin) * sqrt ( (z,z)/(u,u) )
!    (9)  (emax/emin) * sqrt ( (z,z)/(inv(ql)*rhs,inv(ql)*rhs) )
!   (10)                sqrt ( (z,z)/(inv(ql)*rhs,inv(ql)*rhs) )
!
!  Parameters:
!
!        icall      key for initial or secondary call
!                    = 1   initial call
!                    = 2   later call (needed if q is changed)
!        subq       preconditioning routine
!        n          order of system
!        rhs        right hand side
!        ubar       known solution
!        wksp       workspace vector of length n
!        ier        error code
!                    =  0  no error detected
!                    = -7  q is not positive definite
!
!  
!
  dimension rhs(1), ubar(1), wksp(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external subq
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
!
  n = nn
  nteste = ntest
  if (ntest > 6) nteste = ntest - 6
  go to (10,50,20,20,30,40), nteste
!
!  bnorm1: sqrt(b,q(inv)b).
!
 10   call subq (coef,jcoef,wfac,jwfac,n,rhs,wksp)
  sum = vdot (n,rhs,wksp)
  if (sum >= 0.0) go to 15
  ier = -7
  call ershow (ier,'nmcalc')
  return
 15   bnorm1 = amax1 ( sqrt(sum),srelpr )
  return
!
!  bnorm1: sqrt(q(inv)b,q(inv)b).
!
 20   call subq (coef,jcoef,wfac,jwfac,n,rhs,wksp)
  sum = vdot (n,wksp,wksp)
  bnorm1 = amax1 ( sqrt(sum),srelpr )
  return
!
!  bnorm.
!
 30   if (icall == 2) return
  sum = vdot (n,rhs,rhs)
  bnorm = amax1 ( sqrt(sum),srelpr )
  bnorm1 = bnorm
  return
!
!  ubarnm.
!
 40   if (icall == 2) return
  sum = vdot (n,ubar,ubar)
  ubarnm = amax1 ( sqrt(sum),srelpr )
  return
!
!  exit.
!
 50   return
end
subroutine noadp (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! NOADP is a dummy routine to do no adaption.
!
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  return
end
subroutine nspcg (precon,accel,ndimm,mdimm,nn,maxnzz,coef,jcoef,p,ip,u, &
  ubar,rhs,wksp,iwksp,nw,inw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! NSPCG is the driver for the NSPCG package.
!
!  Parameters:
!
!       precon    preconditioning module
!       accel     acceleration module
!       coef      real matrix data array
!       jcoef     integer matrix data array
!       n         input integer.  order of the system (= nn)
!       u         input/output vector.  on input, u contains the
!                 initial guess to the solution.  on output, it
!                 contains the latest estimate to the solution.
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       wksp      real workspace vector of length nw
!       iwksp     integer workspace vector of length inw
!       nw        length of wksp upon input, amount used upon
!                  output
!       inw       length of iwksp upon input, amount used upon
!                  output
!       iparm     integer vector of length 30.  allows user to
!                 specify some integer parameters which affect
!                 the method.
!       rparm     real vector of length 30.  allows user to
!                 specify some real parameters which affect
!                 the method.
!       ier       output integer.  error flag.
!
!  
!
  external  accel, precon
  integer   iparm(30), jcoef(2), p(1), ip(1), iwksp(1)
  dimension coef(1), rhs(1), u(1), ubar(1), rparm(30), wksp(1)
!
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  ier = 0
  ndim = ndimm
  mdim = mdimm
  n = nn
  maxnz = maxnzz
  lenr = nw
  leni = inw
  irmax = 0
  iimax = 0
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,1,ier)
  if (ier < 0) return
  timfac = 0.0
  call pointr (1,wksp,iwksp,ier)
!
!  call preparatory routines.
!
!  remove zeros from jcoef for Purdue data structure.
!
  if (nstore == 1) call adjust (n,ndim,maxnz,jcoef,1)
  call prep (coef,jcoef,wksp(irpnt),iwksp(iipnt),n,nstore,ier)
  if (ier < 0) then
     call ershow (ier,'nspcg')
     go to 20
  end if
!
!  eliminate penalty-method dirichlet points, if requested.
!
  ielim = iparm(24)
  tol = rparm(15)
  if (ielim == 1) call elim (n,jcoef,coef,rhs,wksp,iwksp,tol)
!
!  determine symmetry of matrix.
!
  if (nstore == 1 .and. isymm == 2) call detsym(ndim,maxnz,coef,jcoef,n,isymm)
!
!  scale matrix.
!
  call scale (coef,jcoef,wksp,1,n,u,ubar,rhs,ier)
  if (ier < 0) go to 20
!
!  permute matrix.
!
  call permut (coef,jcoef,p,ip,wksp,iwksp,1,n,u,ubar,rhs,ier)
  if (ier < 0) go to 15
!
!  call iterative routine.
!
  call precon (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!  unpermute matrix.
!
  call permut (coef,jcoef,p,ip,wksp,iwksp,2,n,u,ubar,rhs,ier)
!
!  unscale matrix.
!
 15   call scale (coef,jcoef,wksp,2,n,u,ubar,rhs,ier)
!
!  restore zeros to jcoef for Purdue data structure.
!
 20   if (nstore == 1) call adjust (n,ndim,maxnz,jcoef,2)
  t2 = timer (dummy)
  timtot = t2 - t1
  iparm(18) = ipropa
  iparm(23) = isymm
  rparm(13) = timfac
  rparm(14) = timtot
  call echall (n,iparm,rparm,2,1,ier)
!
  call pointr (2,wksp,iwksp,ier)
  nw = irmax
  inw = iimax
  maxnzz = maxnz
  return
end
subroutine nullpl (coef,jcoef,wk,iwk,n,subql,suba,subqr,u,v)
!
!*******************************************************************************
!
!! NULLPL applies the left preconditioner.
!
  dimension u(1), v(1), coef(1), jcoef(2), wk(1), iwk(1)
  external subql, suba, subqr
!
  call subql (coef,jcoef,wk,iwk,n,u,v)
  return
end
subroutine nullpr (coef,jcoef,wk,iwk,n,subql,suba,subqr,u,v)
!
!*******************************************************************************
!
!! NULLPR applies the right preconditioner.
!
  dimension u(1), v(1), coef(1), jcoef(2), wk(1), iwk(1)
  external subql, suba, subqr
!
  call subqr (coef,jcoef,wk,iwk,n,u,v)
  return
end
subroutine odir (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef, &
  n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! ODIR is the user interface to the ORTHODIR algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call odirw (suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs, &
    wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine odirw (suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
  wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! ODIRW implements ORTHODIR.
!
!
!  The algorithm includes truncation, restarting and 2-sided preconditioning.  
!  the effective value of the z matrix is (inv(ql)*a*inv(qr))**t.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  logical iql, iqr
  external suba, subql, subqr
  dimension iparm(30), rparm(30)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! the following indexing functions are used to access the old
! direction vectors and dot products --
!
  indpt(i) = ipt + mod(i,nv)*n
  indqap(i) = iqapt + mod(i,nv)*n
  inddot(i) = idot + mod(i,nv)
!
! various preliminary calculations.
!
!
  nwusd = 0
  ier = 0
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  write (nout,496)
496   format (' orthodir')
  iacel = 7
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
!
! initialize the stopping test.
!
  call inithv (0)
  zhave  = .true.
  zthave = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 997
!
! memory allocation, etc.
!
  nv = max(1,min(ns1,ns2-1))
  ipt = 1
  iqapt = ipt + nv*n
  idot = iqapt + nv*n
  iz = idot + nv
  izt = iz + n
  if (.not. iqr) izt = iz
  isv = izt + n
  iv1 = isv + n
  iv2 = iv1 + n
!
  if (iql) nwusd = max(nwusd,iv2-1+n)
  if (.not. iql) nwusd = max(nwusd,iv1-1+n)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
!
! perform first-iterate calculations
!
  if (iql) go to 122
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iz))
  call vexopy (n,wk(iz),rhs,wk(iz),2)
  go to 121
 122  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iz))
 121  if (iqr) call subqr (coef,jcoef,wfac,jwfac,n,wk(iz),wk(izt))
!     if (.not. iqr) zdot = vdot (n,wk(iz),wk(iz))
!
!  Begin iteration loop.
!
!
! determine whether or not to stop.
!
 10   call inithv (1)
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,wk(iz),wk(izt),wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
! proceed to calculate the direction vectors.
!
! first, case of no old p vectors.
!
  np = min(mod(in,ns2),ns1)
  if (np /= 0) go to 100
!
  if (is == 0) call vcopy (n,wk(izt),wk(indpt(in)))
  if (is /= 0) call vcopy (n,wk(isv),wk(indpt(in)))
  if (iql) go to 123
  call suba (coef,jcoef,wfac,jwfac,n,wk(indpt(in)),wk(indqap(in)))
  go to 120
 123  call suba (coef,jcoef,wfac,jwfac,n,wk(indpt(in)),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(indqap(in)))
  go to 120
!
! case of at least one old p vector.
! this case is handled in a tricky way, to optimize the workspace.
!
 100  if (iql) go to 124
  call suba (coef,jcoef,wfac,jwfac,n,wk(isv),wk(iv1))
  go to 125
 124  call suba (coef,jcoef,wfac,jwfac,n,wk(isv),wk(iv2))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv2),wk(iv1))
!
 125  top = vdot (n,wk(indqap(in-np)),wk(iv1))
  bet = - top / wk(inddot(in-np))
  call vtriad (n,wk(indpt(in)),wk(isv),bet,wk(indpt(in-np)),1)
  call vtriad (n,wk(indqap(in)),wk(iv1),bet,wk(indqap(in-np)),1)
  ibegin = in - np + 1
  iend = in - 1
  if (ibegin > iend) go to 613
  do 612 i = ibegin,iend
  top = vdot (n,wk(indqap(i)),wk(iv1))
  bet = - top / wk(inddot(i))
  call vtriad (n,wk(indpt(in)),wk(indpt(in)),bet,wk(indpt(i)),1)
 612  call vtriad (n,wk(indqap(in)),wk(indqap(in)),bet,wk(indqap(i)),1)
 613  continue
!
! periodically scale the direction vector, to prevent overflow.
!
 120  continue
  dot = vdot (n,wk(indqap(in)),wk(indqap(in)))
  if (dot<srelpr**2 .or. dot>(1e0/srelpr**2)) then
    call vtriad (n,wk(indpt(in)), xxx,1e0/dot,wk(indpt(in)), 2)
    call vtriad (n,wk(indqap(in)),xxx,1e0/dot,wk(indqap(in)),2)
    dot = 1e0
  end if
!
! at this point, we are finished forming the latest direction vector.
! we proceed to calculate lambda and update the solution and
! the residuals.
!
 129  continue
!     if (abs(dot) < srelpr) go to 998
  wk(inddot(in)) = dot
  top = vdot (n,wk(indqap(in)),wk(iz))
  vlamda = top / dot
! the following commented-out line is unstable.  but it can be fixed.
!     if (.not. iqr) zdot = zdot - 2*vlamda*top + vlamda**2*dot
!
! u --
!
  call vtriad (n,u,u,vlamda,wk(indpt(in)),1)
!
! z --
!
  call vtriad (n,wk(iz),wk(iz),-vlamda,wk(indqap(in)),1)
!
! zt --
!
  call subqr (coef,jcoef,wfac,jwfac,n,wk(indqap(in)),wk(isv))
  if (iqr) call vtriad (n,wk(izt),wk(izt),-vlamda,wk(isv),1)
!
! proceed to next iteration
!
  in = in + 1
  is = is + 1
  if (is == ns2) is = 0
  go to 10
!
!  Finish up.
!
 900  if (.not. halt) go to 996
  if (level >= 1) write (nout,720) in
 720  format (/' orthodir converged in ',i5,' iterations.')
!
 725  continue

  if (idgts >= 0) then
    call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
  end if
!
! pack revised parms into iparm, rparm.
  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
!
! error returns
!
! no convergence.
 996  ier = 1
  call ershow (ier,'odirw')
  zeta = stptst
  go to 725
!
! generic error handler.
 997  call ershow (ier,'odirw')
  go to 735
!
! breakdown.
 998  ier = -15
  call ershow (ier,'odirw')
  go to 725
!
! insufficient real wksp.
 999  ier = -2
  call ershow (ier,'odirw')
  go to 735
end
subroutine omgchg (ssorcp,coef,jcoef,wfac,jwfac,n,p,r)
!
!*******************************************************************************
!
!! OMGCHG changes ALPHAB and BETAB for a new estimate of OMEGA.
!
!
!  Parameters:
!
!         n       order of system (= nn)
!         p       vector from acceleration algorithm
!         r       workspace vector from acceleration algorithm
!
!  
!
  dimension p(1), r(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external ssorcp
!
!
!
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
!
!  update alphab and betab.
!
  call ssorcp (coef,jcoef,wfac,jwfac,n,p,r,pdp,pldup)
  alphab = amin1 (alphab, (pap/pdp) - 1.0)
  betab  = amax1 (betab , pldup/pdp)
  return
end
subroutine omin (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef, &
  n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! OMIN is the user interface to the truncated/restarted ORTHOMIN algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call ominw (suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs, &
    wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine omingw (suba,subql,subqr,precl,precr,coef,jcoef,wfac,jwfac,n, &
  u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! OMINGW is a generalized version of the OMINW routine.
!
!  It allows a more general computational form for the preconditioning.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  logical ipl, ipr
  external suba, subql, subqr, precl, precr
  dimension iparm(30), rparm(30)
  logical ztget, havest, hadest, evest
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! the following indexing functions are used to access the old
! direction vectors and dot products --
!
  indpt(i) = ipt + mod(i,nv)*n
  indqap(i) = iqapt + mod(i,nv)*n
  inddot(i) = idot + mod(i,nv+1)
  indhes(i,j) = ihess + (i-1) + (j-1)*nhess
  inapar(i) = iapar + mod(i,nv)
  indlam(i) = ilam + mod(i,nv+1)
!
! various preliminary calculations.
!
  t1 = timer (dummy)
!
  ipl = iplr == 1 .or. iplr == 3
  ipr = iplr == 2 .or. iplr == 3
!
  iacel = 8
  nwusd = 0
  if (level >= 1) write (nout,497)
497   format (' omin')
!
! initialize the stopping test.
!
  call inithv (0)
  zhave   = .true.
  zthave  = .true.
  nwpstp =  nw
  call pstopg (0,suba,subql,subqr,precl,precr,coef,jcoef,wfac,jwfac,n,u, &
    ubar,rhs,xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 997
  ztget = ztcalp
  zthave = ztget
!
! memory allocation, etc.
!
  numbig = 1000
  methev = 1
  if (iabs(ns3) >= numbig) then
    if (ns3 > 0) ns3 = ns3 - numbig
    if (ns3 < 0) ns3 = ns3 + numbig
    methev = 2
  end if
!
  evest = ns3/=0 .and. (maxadd.or.minadd)
  nhess = 2 + min(itmax,ns2)
  nv = max(1,min(ns1,ns2-1))
  ipt = 1
  iqapt = ipt + nv*n
  idot = iqapt + nv*n
  iapar = idot + (nv+1)
  ihess = iapar + nv
  ilam = ihess + nhess*(nv+2)
  if (.not. evest) ilam = ihess
  iz = ilam + (nv+1)
  izt = iz + n
  if (.not. ipr) izt = iz
  iv1 = izt + n
  iv2 = iv1 + n
  ir = iz
  if (ipl) ir = iv1
!
  nwtmp = iv1 - 1 + n
  if (ipl) nwtmp = iv2 - 1 + n
  nwusd = max(nwusd,nwtmp)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
!
! perform first-iterate calculations
!
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(ir))
  call vexopy (n,wk(ir),rhs,wk(ir),2)
  if (ipl) call precl (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,wk(ir),wk(iz))
  hadest = .false.
!
!  Begin iteration loop.
!
! determine whether or not to stop.
!
 10   if (.not. ztget) go to 710
  if (ipr) call precr (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,wk(iz),wk(izt))

 710  call inithv (1)
  nwpstp = nw - (iv1-1)
  call pstopg (1,suba,subql,subqr,precl,precr,coef,jcoef,wfac,jwfac,n,u, &
    ubar,rhs,xxx,wk(iz),wk(izt),wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
  if (zthave) go to 711
  if (ipr) call precr (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,wk(iz),wk(izt))
!
!  Proceed to calculate the direction vectors.
!
! first, case of no old p vectors.
!
 711  np = min(mod(in,ns2),ns1)
  if (np /= 0) go to 100
!
  call vcopy (n,wk(izt),wk(indpt(in)))
  if (.not. ipl) then
    call suba (coef,jcoef,wfac,jwfac,n,wk(indpt(in)),wk(indqap(in)))
  else
    call suba (coef,jcoef,wfac,jwfac,n,wk(indpt(in)),wk(iv1))
    call precl (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,wk(iv1),wk(indqap(in)))
  end if
  go to 120
!
! case of at least one old p vector.
! this case is handled in a tricky way, to optimize the workspace.
!
 100  if (.not. ipl) then
    call suba (coef,jcoef,wfac,jwfac,n,wk(izt),wk(iv1))
  else
    call suba (coef,jcoef,wfac,jwfac,n,wk(izt),wk(iv2))
    call precl (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,wk(iv2),wk(iv1))
  end if

  top = vdot (n,wk(indqap(in-np)),wk(iv1))
  wk(inapar(in-np)) = top
  bet = - top / wk(inddot(in-np))
  call vtriad (n,wk(indpt(in)),wk(izt),bet,wk(indpt(in-np)),1)
  call vtriad (n,wk(indqap(in)),wk(iv1),bet,wk(indqap(in-np)),1)

  do 612 i = in-np+1, in-1
  top = vdot (n,wk(indqap(i)),wk(iv1))
  wk(inapar(i)) = top
  bet = - top / wk(inddot(i))
  call vtriad (n,wk(indpt(in)), wk(indpt(in)), bet,wk(indpt(i)), 1)
 612  call vtriad (n,wk(indqap(in)),wk(indqap(in)),bet,wk(indqap(i)),1)
!
! at this point, we are finished forming the latest direction vector.
! we proceed to calculate lambda and update the solution and
! the residuals.
!
120   continue
  apap = vdot (n,wk(indqap(in)),wk(indqap(in)))
!     if (abs(apap) < srelpr**2) go to 998
  if (abs(apap) == 0e0) go to 998
  wk(inddot(in)) = apap
  top = vdot (n,wk(indqap(in)),wk(iz))
  vlamda = top / apap
!     if (.not. ipr) zzdot = zzdot - 2*vlamda*top + vlamda**2*apap
!
! u --
  call vtriad (n,u,u,vlamda,wk(indpt(in)),1)
!
! z --
  call vtriad (n,wk(iz),wk(iz),-vlamda,wk(indqap(in)),1)
!
!  Hess matrix update.
!
! there are two schemes here, based on two different ways of projecting
! the iteration matrix.
!
! update Hessenberg matrix: scheme 1
!
  if (.not. evest) go to 955
  wk(indlam(in)) = vlamda
  if (is == 0) call vfill (nhess*(nv+2),wk(ihess),0e0)
  if (methev /= 1) go to 746
!
  do 954 i=in-np,in
  if (i == in) apar = apap
  if (i /= in) apar = wk(inapar(i))
  wk(indhes(i+1+(is-in),in-i+2)) = wk(indhes(i+1+(is-in),in-i+2)) &
    + apar/wk(indlam(in)) / sqrt(wk(inddot(in))*wk(inddot(i)))
  if (is /= 0) &
    wk(indhes(i+1+(is-in),in-i+1)) = wk(indhes(i+1+(is-in),in-i+1)) &
    - apar/wk(indlam(in-1)) / sqrt(wk(inddot(in-1))*wk(inddot(i)))
 954  continue
  iesize = is
  go to 747
!
! update Hessenberg matrix: scheme 2
!
 746  iesize = is + 1
  wk(indhes(is+2,1)) = -1e0 / vlamda
  wk(indhes(is+1,2)) =  1e0 / vlamda
  if (np == 0) go to 749
  do 748 i=in-np,in-1
  id = in - i + 1
  wk(indhes(is+3-id,id  )) = wk(indhes(is+3-id,id  )) &
    - wk(inapar(i))/wk(inddot(i))/wk(indlam(i))
 748  wk(indhes(is+2-id,id+1)) = wk(indhes(is+2-id,id+1)) &
   + wk(inapar(i))/wk(inddot(i))/wk(indlam(i))
 749  continue
!
! estimate eigenvalues.
!
 747  nwhe = nw - (iv1-1)
  call hesest (wk(ihess),nhess,nv+2,iesize,ns3,havest,emaxnw,eminnw, &
    wk(iv1),nwhe,ier)
  nwusd = max (nwusd,iv1-1+nwhe)
  if (ier /= 0) go to 995
  if (.not. havest) go to 955
  if (hadest) go to 956
  if (maxadd) emax = emaxnw
  if (minadd) emin = eminnw
  hadest = .true.
  go to 955
 956  if (maxadd) emax = amax1 (emax,emaxnw)
  if (minadd) emin = amin1 (emin,eminnw)
!
!  Proceed to next iteration.
!
 955  in = in + 1
  is = is + 1
  if (is == ns2) is = 0
  go to 10
!
!  Finish up.
!
 900  if (.not. halt) go to 996
  if (level >= 1) write (nout,720) in
 720  format (/' orthomin converged in ',i5,' iterations.')
!
 725  if (idgts >= 0) then
       call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
  end if
!
! pack revised parms into iparm, rparm.
  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
!  Error returns.
!
! unimplemented option.
 995  ier = -16
  call ershow (ier,'omingw')
  go to 725
!
! no convergence.
 996  ier = 1
  call ershow (ier,'omingw')
  zeta = stptst
  go to 725
!
! generic error handler.
 997  call ershow (ier,'omingw')
  go to 735
!
! breakdown.
 998  ier = -15
  call ershow (ier,'omingw')
  go to 725
!
! insufficient real wksp.
 999  ier = -2
  call ershow (ier,'omingw')
  go to 735
end
subroutine ominw (suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs,wk, &
  nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! OMINW implements the truncated/restarted ORTHOMIN algorithm.
!
!
! eigenvalue estimation is implemented.
! note that this also implements the gcr algorithm.
!
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subql, subqr
  external nullpl, nullpr
  dimension iparm(30), rparm(30)
!
  ier = 0
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) return
!
! pass on to workhorse routine.
!
  call omingw (suba,subql,subqr,nullpl,nullpr,coef,jcoef,wfac,jwfac,n,u, &
    ubar,rhs,wk,nw,iparm,rparm,ier)
  return
end
subroutine ores (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! ORES is the user interface to the ORTHORES algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2), wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call oresw (suba,subql,subqr,coef,jcoef,wksp,iwksp,n,u,ubar,rhs, &
    wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine oresw (suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
  wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! ORESW implements ORTHORES.
!
!
!  The algorithm includes truncation, restarting and 2-sided preconditioning.  
!  the value of z is
! the identity.  the code is optimal in speed and workspace
! requirements, for general a, ql and qr.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subql, subqr
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! the following indexing functions are used to access the old
! direction vectors and dot products --
!
  indu(i) = iu + mod(i,nv)*n
  indz(i) = iz + mod(i,nv)*n
  inddot(i) = idot + mod(i,nv)
!
! various preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 9
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 997
  if (level >= 2) write (nout,496)
496   format (' orthores')
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
!
! initialize the stopping test.
!
  call inithv (0)
  zhave  = .true.
  zthave = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
! memory allocation, etc.
!
! nomenclature -- r  -- residual of the original system.
!                 z  -- inv(ql)*r
!                 zt -- inv(qr)*z
!
  nv = max(1,min(ns1+1,ns2))
  iu = 1
  iz = iu + nv*n
  idot = iz + nv*n
  iv1 = idot + nv
  iv2 = iv1 + n
  nwusd = max(nwusd,iv2-1+n)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
!
! perform first-iterate calculations.
! note -- we will use the vector 'u' to store ztilde. the u vectors
! will be stored in the table.  wk(iv1) will hold r.
!
  call vcopy (n,u,wk(indu(0)))
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(indz(0)))
  call subqr (coef,jcoef,wfac,jwfac,n,wk(indz(0)),u)
  wk(inddot(0)) = vdot (n,wk(indz(0)),wk(indz(0)))
!
!  Begin iteration loop.
!
! determine whether or not to stop --
!
 10   call inithv (1)
  nwpstp = nw - (iv2-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,wk(indu(in)), &
    ubar,rhs,xxx,wk(indz(in)),u,wk(iv2),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv2-1)
  if (level >= 2) call iterm (n,wk(indu(in)))
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
! proceed to calculate the iterates.
!
  np = min(mod(in,ns2)+1,ns1+1)
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  top = vdot (n,wk(indz(in+1-np)),wk(iv2))
  sig = top / wk(inddot(in+1-np))
  call vtriad (n,wk(indz(in+1)),wk(iv2),-sig,wk(indz(in+1-np)),1)
  call vtriad (n,wk(indu(in+1)),u,sig,wk(indu(in+1-np)),1)
  sigsum = sig
  ibegin = in - np + 2
  iend = in
  if (ibegin > iend) go to 613
  do 612 i = ibegin,iend
  top = vdot (n,wk(indz(i)),wk(iv2))
  den = wk(inddot(i))
  if (abs(den) < srelpr) go to 998
  sig = top / den
  call vtriad (n,wk(indz(in+1)),wk(indz(in+1)),-sig,wk(indz(i)),1)
  call vtriad (n,wk(indu(in+1)),wk(indu(in+1)),sig,wk(indu(i)),1)
 612  sigsum = sigsum + sig
 613  continue
  if (abs(sigsum) < srelpr) go to 998
  vlamda = 1.0/sigsum
  call vtriad (n,wk(indz(in+1)),xxx,-vlamda,wk(indz(in+1)),2)
  call vtriad (n,wk(indu(in+1)),xxx,vlamda,wk(indu(in+1)),2)
  wk(inddot(in+1)) = vdot (n,wk(indz(in+1)),wk(indz(in+1)))
!
  call subqr (coef,jcoef,wfac,jwfac,n,wk(indz(in+1)),u)
!
! proceed to next iteration
!
  in = in + 1
  go to 10
!
!  Finish up.
!
 900  call vcopy (n,wk(indu(in)),u)
  if (halt) go to 715
  ier = 1
  call ershow (ier,'oresw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' orthores converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 997  call ershow (ier,'oresw')
  go to 735
!
 998  ier = -15
  call ershow (ier,'oresw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'oresw')
  go to 735
!
end
subroutine out (nn,v,iswt,noutt)
!
!*******************************************************************************
!
!! OUT prints the residual and solution vectors.
!
!
!  Parameters:
!
!          v      vector of length n
!          iswt   labelling information
!          nout   output device number (= noutt)
!
!  
!
  dimension v(nn)
!
     n = nn
     nout = noutt
     if (n <= 0) return
!
     kupper = min (n, 4)
     if (iswt == 1) write (nout,10)
 10      format (//5x,'residual vector')
     if (iswt == 2) write (nout,15)
 15      format (//5x,'solution vector')
     write (nout,20) (i,i=1,kupper)
 20      format (10x,4i15)
     write (nout,25)
 25      format (10x,65('-') /)
!
     do 35 j = 1,n,4
        kupper = min (j+3,n)
        jm1 = j - 1
        write (nout,30) jm1,(v(k),k=j,kupper)
 30         format (4x,i5,'+  ',4e15.5)
 35      continue
!
     return
end
subroutine parsi
!
!*******************************************************************************
!
!! PARSI computes the iteration parameters.
!
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  rhol = rho
  if (is - 1) 10,15,20
 10   rho = 1.0
  go to 25
 15   rho = 1.0/(1.0 - sigma*sigma/2.0)
  go to 25
 20   rho = 1.0/(1.0 - sigma*sigma*rhol/4.0)
!
!  compute alpha, beta.
!
 25   alpha = rho*gamma
  beta = rhol*(rho - 1.0)/rho
  return
end
subroutine pbneu (suba,dsolve,coef,jcoef,wfac,jwfac,nd,wksp,nn,r,z)
!
!*******************************************************************************
!
!! PBNEU computes a block Neumann polynomial approximation to inverse(A).
!
!
!     if a = d - b, where d is a dense banded matrix
!     then the output vector is --
!
!             z = (i + p + p**2 +. + p**nd)*inv(d) * r
!
!     where   p = inv(d)*b .
!
!  Parameters:
!
!         suba    matrix-vector multiplication routine
!         dsolve  routine for computing inv(d)*vector
!         nd      the degree of the polynomial desired
!         wksp    workspace of length 2*n
!         n       order of system (= nn)
!         r       residual
!         z       output vector
!
!  
!
  external  suba, dsolve
  dimension r(1), z(1), wksp(1), coef(1), jcoef(2), wfac(1), jwfac(1)
!
  n = nn
  np1 = n + 1
  call dsolve (coef,jcoef,wfac,jwfac,n,r,z)
  if (nd <= 0) return
!
  do 20 k = 1,nd
     call suba (coef,jcoef,wfac,jwfac,n,z,wksp)
     do 10 i = 1,n
 10      wksp(i) = r(i) - wksp(i)
     call dsolve (coef,jcoef,wfac,jwfac,n,wksp,wksp(np1))
     do 15 i = 1,n
 15      z(i) = z(i) + wksp(i+n)
 20   continue
  return
end
subroutine pbpii (suba,dsolve,coef,jcoef,wfac,jwfac,ainf,alpha,beta,nd, &
  wksp,nn,r,z)
!
!*******************************************************************************
!
!! PBPII computes a block least squares polynomial approximation to inverse(A).  
!
!  the output vector is --
!
!         z = inv(d)*p  (a*inv(d)) * r
!                     np
!
!  Parameters:
!
!         suba    matrix-vector multiplication routine
!         dsolve  routine to compute inv(d)*vector
!         ainf    the infinity norm of matrix inv(d)*a
!         alpha,  the least squares weighting factors
!          beta
!         nd      the degree of the polynomial desired
!         wksp    workspace of length 2*n
!         n       order of system (= nn)
!         r       residual
!         z       output vector
!
!  
!
  external  suba, dsolve
  dimension r(1), z(1), wksp(1), coef(1), jcoef(2), wfac(1), jwfac(1)
!
!
  n = nn
  np1 = n + 1
  al = alpha
  be = beta
!
  c1 = ((al+be+2.0)*(al+be+3.0))/(ainf*(al+2.0)*(al+be+2.0))
  call dsolve (coef,jcoef,wfac,jwfac,n,r,z)
  do 10 i = 1,n
 10   z(i) = c1*z(i)
  if (nd <= 0) return
!
  do 15 i = 1,n
 15   wksp(i) = r(i)
  do 35 k = 1,nd
     fk = float (k)
     c1 = ((2.0*fk+al+be+2.0)*(2.0*fk+al+be+3.0)) &
       / (ainf*(fk+al+2.0)*(fk+al+be+2.0))
     c2 = (fk*(fk+be)*(2.0*fk+al+be))/ &
       ((fk+al+1.0)*(fk+al+be+1.0)*(2.0*fk+al+be+2.0))
     call suba (coef,jcoef,wfac,jwfac,n,z,wksp(np1))
     do 20 i = 1,n
 20      wksp(n+i) = r(i) - wksp(n+i)
     do 25 i = 1,n
 25      wksp(i) = wksp(i+n) + c2*wksp(i)
     call dsolve (coef,jcoef,wfac,jwfac,n,wksp,wksp(np1))
     do 30 i = 1,n
 30      z(i) = z(i) + c1*wksp(n+i)
 35   continue
  return
end
subroutine pbs (n,t1,t2,x)
!
!*******************************************************************************
!
!! PBS does a penta-diagonal back substitution.
!
!
!  This has the form (i+t1+t2)*x = y
!     where t1 and t2 are the first and second super diagonals.
!
!  Parameters:
!
!          n      order of the system
!          t1     vector of length n-1 containing the first super-
!                  diagonal elements
!          t2     vector of length n-2 containing the second super-
!                  diagonal elements
!          x      on input, x contains y
!                 on output, x contains the solution to
!                  (i + t1 + t2)*x = y
!
!  
!
  dimension t1(1), t2(1), x(1)
!
  x(n-1) = x(n-1) - t1(n-1)*x(n)
  do 10 i = n-2,1,-1
 10   x(i) = x(i) - t1(i)*x(i+1) - t2(i)*x(i+2)
  return
end
subroutine pbsm (nn,nsize,t1,t2,x)
!
!*******************************************************************************
!
!! PBSM does a penta-diagonal back substitution.
!
!  This has the form (i+t1+t2)*x = y
!     where t1 and t2 are superdiagonals of a system composed of
!     independent subsystems of size nsize.
!
!  Parameters:
!
!          n      order of system
!          nsize  order of the individual subsystems
!          t1     linear array of length n-1 containing the first
!                  super-diagonal elements of the factorizations
!          t2     linear array of length n-2 containing the second
!                  super-diagonal elements of the factorizations
!          x      on input, x contains y
!                 the solution to (i + t1 + t2)*x = y
!
!  
!
  dimension t1(nsize,1), t2(nsize,1), x(nsize,1)
!
  n = nn
  nsys = n/nsize
  do 10 j = 1,nsys
 10   x(nsize-1,j) = x(nsize-1,j) - t1(nsize-1,j)*x(nsize,j)
  do 20 i = nsize-2,1,-1
     do 15 j = 1,nsys
 15      x(i,j) = x(i,j) - t1(i,j)*x(i+1,j) - t2(i,j)*x(i+2,j)
 20   continue
  return
end
subroutine permas (isym,nn,nzz,ia,ja,a,wksp,p)
!
!*******************************************************************************
!
!! PERMAS permutes the rows and columns of a sparse matrix.
!
!
!  Parameters:
!
!          isym      switch for symmetric storage
!                     = 0   matrix is symmetric
!                     = 1   matrix is nonsymmetric
!          n         size of system
!          nz        length of ia, ja, and a vectors
!          ia        vector of i values
!          ja        vector of j values
!          a         vector of matrix coefficients
!          wksp      workspace vector of length n
!          p         permutation vector
!
!  it is assumed that the i-th entry of the permutation vector
!     p indicates the row the i-th row gets mapped into.  (i.e.
!     if ( p(i) = j ) row i gets mapped into row j)
!
!  
!
  dimension a(1), wksp(1)
  integer   ia(1), ja(1), p(1)
!
  n = nn
  nz = nzz
!
!  explicit gathers.
!
  call vgathi (nz,p,ia,ia)
  call vgathi (nz,p,ja,ja)
  do 5 i = 1,n
 5    wksp(i) = a(i)
  call vscatr (n,wksp,p,a)
  do 10 i = 1,n
     ia(i) = i
     ja(i) = i
 10   continue
!
!  convert to upper triangular elements for symmetric storage
!
  if (isym == 1) return
  np1 = n + 1
  do 15 i = np1,nz
     if (ia(i) <= ja(i)) go to 15
     idum = ia(i)
     ia(i) = ja(i)
     ja(i) = idum
 15   continue
  return
end
subroutine permat (ndim,maxnz,coef,jcoef,wksp,iwksp,nn,p)
!
!*******************************************************************************
!
!! PERMAT permutes the rows and columns of a Purdue sparse matrix.
!
!
!  Parameters:
!
!          ndim      row dimension of coef array in defining routine
!          maxnz     number of columns in coef and jcoef arrays
!          coef      array of matrix coefficients
!          jcoef     array of matrix columns numbers
!          wksp      workspace array of length n
!          iwksp     integer workspace array of length n
!          n         order of system (= nn)
!          p         permutation vector
!
!  it is assumed that the i-th entry of the permutation vector
!     p indicates the row the i-th row gets mapped into.  (i.e.
!     if ( p(i) = j ) row i gets mapped into row j)
!
!  
!
  dimension coef(ndim,1), wksp(1)
  integer   jcoef(ndim,1), iwksp(1), p(1)
!
  n = nn
  if (n <= 0) return
  do 20 j = 1,maxnz
     do 10 i = 1,n
        wksp(i) = coef(i,j)
        iwksp(i) = jcoef(i,j)
 10      continue
     call vscatr (n,wksp,p,coef(1,j))
     call vscati (n,iwksp,p,jcoef(1,j))
     call vgathi (n,p,jcoef(1,j),jcoef(1,j))
 20   continue
  return
end
subroutine permd (coef,jcoef,p,ip,wksp,iwksp,icall,nn,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! PERMD permutes the matrix, U, UBAR, and RHS. (diagonal format)
!
!
!  Parameters:
!
!       icall     key to indicate whether permuting (icall=1) or
!                  unpermuting (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -2   insufficient real space to permute system
!                  =  -3   insufficient integer space to permute
!                          system
!
!  
!
  integer jcoef(2), p(1), ip(1), iwksp(1)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
!
  n = nn
  if (icall == 2) go to 65
!
!  permute system.
!
!  check for sufficient storage to permute matrix
!
  call needw ('permd',0,irpnt,n,ier)
  if (ier < 0) return
  nc = iipnt
  call needw ('permd',1,nc,n,ier)
  if (ier < 0) return
  call pgen (n,p,ip,iwksp(nc),ncolor)
  ipt = nc + ncolor
  ncmax = 0
  do 16 i = nc,ipt-1
     if (ncmax < iwksp(i)) ncmax = iwksp(i)
 16   continue
  call needw ('permd',1,ipt,ncolor+1,ier)
  if (ier < 0) return
  call iptgen (ncolor,iwksp(ipt),iwksp(nc))
  maxnew = ipt + ncolor + 1
  jcnew = maxnew + ncolor
  lbhb = jcnew + ncolor*mdim
  call needw ('permd',1,maxnew,ncolor+ncolor*mdim+n,ier)
  if (ier < 0) return
  isym = nstore - 2
  call pmdg (ndim,mdim,n,maxnz,jcoef,coef,ncolor,iwksp(nc),p,ip,maxd, &
    iwksp(maxnew),iwksp(jcnew),wksp(irpnt),iwksp(lbhb),isym,ier)
  if (ier < 0) then
     call ershow (ier,'permd')
     return
  end if
  lbhb = jcnew + ncolor*maxd
  iblock = lbhb + ncolor
  call move4 (ndim,n,iwksp(maxnew),iwksp(jcnew),coef,ncolor,iwksp(nc), &
    wksp(irpnt),iwksp(lbhb))
  call needw ('permd',1,lbhb,ncolor+3*ncolor*(maxd+1),ier)
  if (ier < 0) return
  call define (ndim,iwksp(maxnew),iwksp(jcnew),coef,ncolor,iwksp(nc), &
    iwksp(iblock),iwksp(lbhb))
  lbhbm = iwksp(lbhb)
  do 45 j = 1,ncolor-1
     if (lbhbm < iwksp(lbhb+j)) lbhbm = iwksp(lbhb+j)
 45   continue
  is1 = iblock + 3*ncolor*lbhbm
  is2 = is1 + ncolor
  call needw ('permd',1,is1,2*ncolor,ier)
  if (ier < 0) return
  call prbblk (ncolor,ncolor,iwksp(iblock),iwksp(lbhb),iwksp(is1), &
    iwksp(is2),propa)
  if (propa) ipropa = 1
  if (.not. propa) ipropa = 0
  iipnt = is1
  call pervec (n,p,rhs,wksp(irpnt))
  call pervec (n,p,u,wksp(irpnt))
  if (ntest == 6) call pervec (n,p,ubar,wksp(irpnt))
  return
!
!  unpermute system.
!
 65   call needw ('permd',1,iipnt,2*n,ier)
  if (ier < 0) return
  isym = nstore - 2
  call unpmdg (ndim,n,maxnz,jcoef,coef,ncolor,iwksp(nc),p,ip,maxd, &
    iwksp(maxnew),iwksp(jcnew),wksp(irpnt),iwksp(iipnt),isym)
  call pervec (n,ip,rhs,wksp(irpnt))
  call pervec (n,ip,u,wksp(irpnt))
  if (ntest == 6) call pervec (n,ip,ubar,wksp(irpnt))
  return
end
subroutine permp (coef,jcoef,p,ip,wksp,iwksp,icall,nn,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! PERMP permutes the matrix, U, UBAR, and RHS.  (Purdue format)
!
!
!  Parameters:
!
!       icall     key to indicate whether permuting (icall=1) or
!                  unpermuting (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -2   insufficient real space to permute system
!                  =  -3   insufficient integer space to permute
!                          system
!
!  
!
  integer jcoef(2), p(1), ip(1), iwksp(1)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
!
  n = nn
  if (icall == 2) go to 40
!
!  permute system.
!
!  check for sufficient storage to permute matrix
!
  call needw ('permp',0,irpnt,n,ier)
  if (ier < 0) return
  nc = iipnt
  call needw ('permp',1,nc,n,ier)
  if (ier < 0) return
  call pgen (n,p,ip,iwksp(nc),ncolor)
  ipt = nc + ncolor
  ncmax = 0
  do 20 i = nc,ipt-1
     if (ncmax < iwksp(i)) ncmax = iwksp(i)
 20   continue
  call needw ('permp',1,ipt,ncolor+1,ier)
  if (ier < 0) return
  call iptgen (ncolor,iwksp(ipt),iwksp(nc))
  iipnt = iipnt + 2*ncolor + 1
  call needw ('permp',1,iipnt,n,ier)
  if (ier < 0) return
  call permat (ndim,maxnz,coef,jcoef,wksp(irpnt),iwksp(iipnt),n,p)
  call needw ('permp',1,iipnt,2*ncolor,ier)
  if (ier < 0) return
  ndt = iipnt
  ndb = iipnt + ncolor
  call move3 (ndim,mdim,n,maxnz,jcoef,coef,iwksp(ndt),iwksp(ndb),ncolor, &
    iwksp(nc),ier)
  iipnt = iipnt + 2*ncolor
  if (ier < 0) then
     call ershow (ier,'permp')
     return
  end if
  call pervec (n,p,rhs,wksp(irpnt))
  call pervec (n,p,u,wksp(irpnt))
  if (ntest == 6) call pervec (n,p,ubar,wksp(irpnt))
  return
!
!  unpermute system.
!
 40   call needw ('permp',0,irpnt,n,ier)
  if (ier < 0) return
  call needw ('permp',1,iipnt,n,ier)
  if (ier < 0) return
  call permat (ndim,maxnz,coef,jcoef,wksp(irpnt),iwksp(iipnt),n,ip)
  call pervec (n,ip,rhs,wksp(irpnt))
  call pervec (n,ip,u,wksp(irpnt))
  if (ntest == 6) call pervec (n,ip,ubar,wksp(irpnt))
  return
end
subroutine perms (coef,jcoef,p,ip,wksp,iwksp,icall,nn,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! PERMS permutes the matrix, U, UBAR, and RHS. (sparse format)
!
!  Parameters:
!
!       icall     key to indicate whether permuting (icall=1) or
!                  unpermuting (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -2   insufficient real space to permute system
!                  =  -3   insufficient integer space to permute
!                          system
!
!  
!
  integer jcoef(2), p(1), ip(1), iwksp(1)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  n = nn
  isym = 0
  if (nstore == 5) isym = 1
  if (icall == 2) go to 10
!
!  permute system.
!
!  check for sufficient storage to permute matrix
!
  call needw ('perms',0,irpnt,n,ier)
  if (ier < 0) return
  call needw ('perms',1,iipnt,n,ier)
  if (ier < 0) return
  call pgen (n,p,ip,iwksp(iipnt),ncolor)
  call permas (isym,n,maxnz,jcoef,jcoef(ndim+1),coef,wksp(irpnt),p)
  call pervec (n,p,rhs,wksp(irpnt))
  call pervec (n,p,u,wksp(irpnt))
  if (ntest == 6) call pervec (n,p,ubar,wksp(irpnt))
  return
!
!  unpermute system.
!
 10   call needw ('perms',0,irpnt,n,ier)
  if (ier < 0) return
  call needw ('perms',1,iipnt,n,ier)
  if (ier < 0) return
  call permas (isym,n,maxnz,jcoef,jcoef(ndim+1),coef,wksp(irpnt),ip)
  call pervec (n,ip,rhs,wksp(irpnt))
  call pervec (n,ip,u,wksp(irpnt))
  if (ntest == 6) call pervec (n,ip,ubar,wksp(irpnt))
  return
end
subroutine permut (coef,jcoef,p,ip,wksp,iwksp,icall,n,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! PERMUT permutes the matrix, U, UBAR, and RHS.
!
!  Parameters:
!
!       icall     key to indicate whether permuting (icall=1) or
!                  unpermuting (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -2   insufficient real space to permute system
!                  =  -3   insufficient integer space to permute
!                          system
!
!  
!
  integer jcoef(2), p(1), ip(1), iwksp(1)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!  data common blocks
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  if (iperm /= 1) return
  go to (5,10,10,15,15), nstore
 5    call permp (coef,jcoef,p,ip,wksp,iwksp,icall,n,u,ubar,rhs,ier)
  return
 10   call permd (coef,jcoef,p,ip,wksp,iwksp,icall,n,u,ubar,rhs,ier)
  return
 15   call perms (coef,jcoef,p,ip,wksp,iwksp,icall,n,u,ubar,rhs,ier)
  return
end
subroutine perror2 (suba,coef,jcoef,wfac,jwfac,nn,u,rhs,wksp,digtt1, &
  digtt2,idgtts)
!
!*******************************************************************************
!
!! PERROR2 computes the residual, R = RHS - A*U.  
!
!  the user
!     also has the option of printing the residual and/or the
!     unknown vector depending on idgts.
!
!  Parameters:
!
!          suba   matrix-vector multiplication routine
!          n      dimension of matrix (= nn)
!          u      latest estimate of solution
!          rhs    right hand side of matrix problem
!          wksp   workspace vector of length n
!          digit1 output - measure of accuracy of stopping test (= digtt1
!          digit2 output - measure of accuracy of solution (= digtt2)
!          idgts   parameter controlling level of output (= idgtts)
!                    if idgts < 1 or idgts > 4, then no output.
!                            = 1, then number of digits is printed, pro-
!                                 vided level >= 1
!                            = 2, then solution vector is printed, pro-
!                                 vided level >= 1
!                            = 3, then residual vector is printed, pro-
!                                 vided level >= 1
!                            = 4, then both vectors are printed, pro-
!                                 vided level >= 1
!
!  
!
  external suba
  dimension rhs(1), u(1), wksp(1), coef(1), jcoef(2),wfac(1), jwfac(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
  idgts = idgtts
  digit1 = 0.0
  digit2 = 0.0
!
  digit1 = -alog10 (abs (srelpr))
  if (stptst > 0.0) digit1 = -alog10 (abs (stptst))
     call suba (coef,jcoef,wfac,jwfac,n,u,wksp)
     call vexopy (n,wksp,rhs,wksp,2)
     rnrm = sqrt ( vdot (n,wksp,wksp) )
     sum = vdot (n,rhs,rhs)
     bnorm = amax1 ( sqrt(sum),srelpr )
     temp = rnrm/bnorm
     if (temp == 0.0) go to 10
     digit2 = -alog10 (abs (temp))
     go to 15
!
 10   digit2 = -alog10 (abs (srelpr))
!
 15      if ((idgts < 1) .or. (level <= 0)) go to 25
     write (nout,20) digit1,digit2
 20      format (/10x,'approx. no. of digits in stopping test =', &
               f5.1,2x,'(digit1)' &
               /10x,'approx. no. of digits in ratio test    =', &
               f5.1,2x,'(digit2)')

     if (idgts <= 1 .or. idgts > 4) go to 25
     if (idgts >= 3) call out (n,wksp,1,nout)
     if (idgts /= 3) call out (n,u,2,nout)
!
 25   continue
  digtt1 = digit1
  digtt2 = digit2
  return
end
subroutine pervec (nn,p,v,wksp)
!
!*******************************************************************************
!
!! PERVEC permutes a vector as dictated by the permutation vector.  
!
!  if p(i) = j, then v(j) gets v(i).
!
!  Parameters:
!
!          n       length of vectors p, v, and wksp  (= nn)
!          p       integer permutation vector
!          v       vector to be permuted
!          wksp    workspace vector of length n
!
!  
!
  integer   p(1)
  dimension v(1), wksp(1)
!
  n = nn
  if (n <= 0) return
  do 10 i = 1,n
     wksp(i) = v(i)
 10   continue
  call vscatr (n,wksp,p,v)
  return
end
subroutine pfac (nn,d,t1,t2)
!
!*******************************************************************************
!
!! PFAC computes a factorization of a single symmetric pentadiagonal matrix.
!
!
!  The matrix is contained in d, t1, and t2 and the factored matrix replaces it.
!
!  Parameters:
!
!          n      order of system (= nn)
!          d      vector of length n containing the diagonal
!                  elements of the matrix
!          t1     vector of length n-1 containing the first
!                  super-diagonal elements of the matrix
!          t2     vector of length n-2 containing the second
!                  super-diagonal elements of the matrix
!
!  
!
  dimension d(1), t1(1), t2(1)
!
  n = nn
  do 10 i = 1,n-2
     dii = 1.0/d(i)
     d(i+1) = d(i+1) - t1(i)*t1(i)*dii
     d(i+2) = d(i+2) - t2(i)*t2(i)*dii
     t1(i+1) = t1(i+1) - t1(i)*t2(i)*dii
 10   continue
  d(n) = d(n) - t1(n-1)*t1(n-1)/d(n-1)
  do 15 i = 1,n
 15   d(i) = 1.0/d(i)
  do 20 i = 1,n-1
 20   t1(i) = d(i)*t1(i)
  do 25 i = 1,n-2
 25   t2(i) = d(i)*t2(i)
  return
end
subroutine pfacm (nn,nsize,d,t1,t2)
!
!*******************************************************************************
!
!! PFACM factors multiple independent symmetric pentadiagonal matrices.
!
!
!  The matrices are contained in d, t1, and t2.
!
!  Parameters:
!
!          n      order of global system (= nn)
!          nsize  size of the individual subsystems
!          d      linear array of length n containing the
!                  diagonal elements of the systems
!          t1     linear array of length n-1 containing the
!                  first super-diagonal elements of the systems
!          t2     linear array of length n-2 containing the
!                  second super-diagonal elements of the systems
!
!  
!
  dimension d(nsize,1), t1(nsize,1), t2(nsize,1)
!
  n = nn
  nsys = n/nsize
  do 15 i = 1,nsize-2
     do 10 j = 1,nsys
        d(i+1,j) = d(i+1,j) - (t1(i,j)**2)/d(i,j)
        d(i+2,j) = d(i+2,j) - (t2(i,j)**2)/d(i,j)
        t1(i+1,j) = t1(i+1,j) - t1(i,j)*t2(i,j)/d(i,j)
 10      continue
 15   continue
  do 20 j = 1,nsys
 20   d(nsize,j) = d(nsize,j) - (t1(nsize-1,j)**2)/d(nsize-1,j)
  call vinv (n,d)
  call vexopy (n-1,t1,d,t1,3)
  call vexopy (n-2,t2,d,t2,3)
  return
end
subroutine pfacn (nn,d,t1,t2,b1,b2)
!
!*******************************************************************************
!
!! PFACN factors a nonsymmetric pentadiagonal matrix.
!
!
!  The matrix is contained in d,t1,t2,b1, and b2 and the factor replaces it.
!
!  Parameters:
!
!          n      order of system (= nn)
!          d      vector of length n containing the diagonal
!                  elements of the matrix
!          t1     vector of length n-1 containing the first
!                  super-diagonal elements of the matrix
!          t2     vector of length n-2 containing the second
!                  super-diagonal elements of the matrix
!          b1     vector of length n-1 containing the first
!                  sub-diagonal elements of the matrix
!          b2     vector of length n-2 containing the second
!                  sub-diagonal elements of the matrix
!
!  
!
  dimension d(1), t1(1), t2(1), b1(1), b2(1)
!
  n = nn
  do 10 i = 1,n-2
     dii = 1.0/d(i)
     d(i+1) = d(i+1) - b1(i)*t1(i)*dii
     d(i+2) = d(i+2) - b2(i)*t2(i)*dii
     t1(i+1) = t1(i+1) - b1(i)*t2(i)*dii
     b1(i+1) = b1(i+1) - b2(i)*t1(i)*dii
 10   continue
  d(n) = d(n) - b1(n-1)*t1(n-1)/d(n-1)
  do 15 i = 1,n
 15   d(i) = 1.0/d(i)
  do 20 i = 1,n-1
     t1(i) = d(i)*t1(i)
     b1(i) = d(i)*b1(i)
 20   continue
  do 25 i = 1,n-2
     t2(i) = d(i)*t2(i)
     b2(i) = d(i)*b2(i)
 25   continue
  return
end
subroutine pfacnm (nn,nsize,d,t1,t2,b1,b2)
!
!*******************************************************************************
!
!! PFACNM factors multiple independent nonsymmetric pentadiagonal matrices.
!
!
!  The matrices are contained in
!     d,t1,t2,b1, and b2.
!
!  Parameters:
!
!          n      order of global system (= nn)
!          nsize  order of single subsystem
!          d      linear array of length n containing the
!                  diagonal elements of the systems
!          t1     linear array of length n-1 containing the first
!                  super-diagonal elements of the systems
!          t2     linear array of length n-2 containing the second
!                  super-diagonal elements of the systems
!          b1     linear array of length n-1 containing the first
!                  sub-diagonal elements of the systems
!          b2     linear array of length n-2 containing the second
!                  sub-diagonal elements of the systems
!
!  
!
  dimension d(nsize,1), t1(nsize,1), b1(nsize,1), t2(nsize,1), b2(nsize,1)
!
  n = nn
  nsys = n/nsize
  do 15 i = 1,nsize-2
     do 10 j = 1,nsys
        d(i+1,j) = d(i+1,j) - b1(i,j)*t1(i,j)/d(i,j)
        d(i+2,j) = d(i+2,j) - b2(i,j)*t2(i,j)/d(i,j)
        t1(i+1,j) = t1(i+1,j) - b1(i,j)*t2(i,j)/d(i,j)
        b1(i+1,j) = b1(i+1,j) - b2(i,j)*t1(i,j)/d(i,j)
 10      continue
 15   continue
  do 20 j = 1,nsys
 20   d(nsize,j) = d(nsize,j) - b1(nsize-1,j)*t1(nsize-1,j)/ d(nsize-1,j)
  call vinv (n,d)
  call vexopy (n-1,t1,d,t1,3)
  call vexopy (n-2,t2,d,t2,3)
  call vexopy (n-1,b1,d,b1,3)
  call vexopy (n-2,b2,d,b2,3)
  return
end
subroutine pfact1 (coef,jcoef,wksp,iwksp,nn,methh,ier)
!
!*******************************************************************************
!
!! PFACT1 computes a point incomplete factorization.
!
!
!  Parameters:
!
!       n       order of system
!       meth    method of factorization
!                = 1   ic   (unmodified)
!                = 2   mic  (modified)
!       nfactr  amount of real workspace needed for factorization
!       nfacti  amount of integer workspace needed for factorization
!       ier     error flag
!
!  
!
!
  integer jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
  n = nn
  meth = methh
!
!  if requested, find out if matrix has property a.
!
  if (ipropa == 0) propa = .false.
  if (ipropa == 1) propa = .true.
  if (lvfill > 0) propa = .false.
  if (lvfill > 0) go to 55
  if (ipropa /= 2) go to 15
  call needw ('pfact1',1,iipnt,2*n,ier)
  if (ier < 0) return
  call prbndx (n,ndim,maxnz,jcoef,coef,iwksp(iipnt),iwksp(iipnt+n),propa,1)
  if (propa) ipropa = 1
  if (.not. propa) ipropa = 0
!
 15   if (.not. propa) go to 35
!
!  propa = .true.
!
  ifactr = irpnt
  nfactr = n
  nfacti = 0
  call needw ('pfact1',0,irpnt,nfactr+n,ier)
  if (ier < 0) return
  call vcopy (n,coef,wksp(ifactr))
  irpnt = irpnt + nfactr
  ip1 = ndim + 1
  ip2 = ndim*(maxt + 1) + 1
  if (isymm == 0) call icfp (ndim,ndim,n,maxt,jcoef(ip1),wksp(ifactr), &
    coef(ip1),meth,1,omega,wksp(irpnt),iflag)
  if (isymm /= 0) call icfnp (ndim,ndim,n,maxt,maxb,jcoef(ip1), &
    jcoef(ip2),wksp(ifactr),coef(ip1),coef(ip2),meth,1,omega,iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfact1')
  return
!
!  propa = .false., lvfill = 0.
!
 35   ifactr = irpnt
  jmax = maxt + 1
  if (isymm /= 0) jmax = 1 + maxt + maxb
  nfactr = n*jmax
  nfacti = 0
  call needw ('pfact1',0,irpnt,nfactr+n,ier)
  if (ier < 0) return
  call vfill (nfactr,wksp(ifactr),0.0)
  do 45 j = 1,jmax
     ip1 = ndim*(j - 1) + 1
     ip2 = n*(j - 1) + ifactr
     call vcopy (n,coef(ip1),wksp(ip2))
 45   continue
  irpnt = irpnt + nfactr
  ip1 = ndim + 1
  ip2 = ndim*(maxt + 1) + 1
  ip3 = ifactr + n
  ip4 = n*(maxt + 1) + ifactr
  if (isymm == 0) call icfp (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr),wksp(ip3), &
    meth,0,omega,wksp(irpnt),iflag)
  if (isymm /= 0) call icfnp (n,ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2), &
    wksp(ifactr),wksp(ip3),wksp(ip4),meth,0,omega,iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfact1')
  return
!
!  propa = .false., lvfill > 0
!
 55   len = n*(maxt + 1)
  if (isymm /= 0) len = n*(1 + maxt + maxb)
  call needw ('pfact1',1,iipnt,len,ier)
  if (ier < 0) return
  call needw ('pfact1',0,irpnt,len,ier)
  if (ier < 0) return
  jmax = maxt + 1
  if (isymm /= 0) jmax = 1 + maxt + maxb
  do 70 j = 1,jmax
     ipt1 = (j - 1)*ndim + 1
     ipt2 = (j - 1)*n + iipnt
     call vicopy (n,jcoef(ipt1),iwksp(ipt2))
     ipt2 = (j - 1)*n + irpnt
     call vcopy (n,coef(ipt1),wksp(ipt2))
 70   continue
  mw1 = (leni - (iipnt + n) + 1)/n
  mw2 = (lenr - (irpnt + n) + 1)/n
  mwidth = min (mw1,mw2)
  maxc = maxt + maxb
  do 75 i = 1,lvfill
     if (isymm == 0) call fillsp (n,n,maxt,iwksp(iipnt+n), &
       wksp(irpnt+n),mwidth,ier)
     if (isymm /= 0) call fillnp (n,n,maxc,iwksp(iipnt+n), &
       wksp(irpnt+n),mwidth,ier)
     if (ier < 0) then
        call ershow (ier,'pfact1')
        return
     end if
 75   continue
  maxcp1 = maxc + 1
  if (isymm /= 0) call move1 (n,mwidth+1,n,maxcp1, &
    iwksp(iipnt),wksp(irpnt),maxt,maxb,ier)
  if (ier < 0) then
     call ershow (ier,'pfact1')
     return
  end if
  if (isymm == 0) nfactr = n*(maxt + 1)
  if (isymm /= 0) nfactr = n*(maxt + maxb + 1)
  nfacti = nfactr
  call needw ('pfact1',0,irpnt,nfactr+n,ier)
  if (ier < 0) return
  call needw ('pfact1',1,iipnt,nfacti,ier)
  if (ier < 0) return
!
  ifactr = irpnt
  ifacti = iipnt
  irpnt = irpnt + nfactr
  iipnt = iipnt + nfacti
  ip1 = ifacti + n
  ip2 = ifacti + n*(maxt + 1)
  ip3 = ifactr + n
  ip4 = ifactr + n*(maxt + 1)
  if (isymm == 0) call icfp (n,n,n,maxt,iwksp(ip1),wksp(ifactr),wksp(ip3), &
    meth,0,omega,wksp(irpnt),iflag)
  if (isymm /= 0) call icfnp (n,n,n,maxt,maxb,iwksp(ip1),iwksp(ip2), &
    wksp(ifactr),wksp(ip3),wksp(ip4),meth,0,omega,iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfact1')
  return
end
subroutine pfact2 (coef,jcoef,wksp,iwksp,nn,methh,ier)
!
!*******************************************************************************
!
!! PFACT2 computes a point incomplete factorization.
!
!  Parameters:
!
!       n       order of system
!       meth    method of factorization
!                = 1   ic   (unmodified)
!                = 2   mic  (modified)
!       nfactr  amount of real workspace needed for factorization
!       ier     error flag
!
!  
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  integer jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
!
  n = nn
  meth = methh
!
!  if requested, find out if matrix has property a.
!
  if (ipropa == 0) propa = .false.
  if (ipropa == 1) propa = .true.
  if (lvfill > 0) propa = .false.
  if (lvfill > 0) go to 20
  if (ipropa /= 2) go to 15
  call needw ('pfact2',1,iipnt,2*n,ier)
  if (ier < 0) return
  call prbndx (n,ndim,maxnz,jcoef,coef,iwksp(iipnt),iwksp(iipnt+n),propa,2)
  if (propa) ipropa = 1
  if (.not. propa) ipropa = 0
!
 15   if (.not. propa) go to 20
!
!  propa = .true.
!
  maxt = maxnz - 1
  maxb = 0
  ifactr = irpnt
  nfactr = n
  nfacti = 0
  call needw ('pfact2',0,irpnt,nfactr+n,ier)
  if (ier < 0) return
  call rowise (maxnz,jcoef,irwise)
  call needw ('pfact2',1,iipnt,maxnz+maxt**2,ier)
  if (ier < 0) return
  call vfill (n,wksp(ifactr),0.0)
  call vcopy (n,coef,wksp(ifactr))
  irpnt = irpnt + nfactr
  if (ifctv == 0) call icf (ndim,n,maxt,jcoef(2),wksp(ifactr),coef(ndim+1), &
    meth,1,omega,wksp(irpnt),iwksp(iipnt),iflag)
  if (ifctv == 1) call icfv(ndim,n,maxt,jcoef(2),wksp(ifactr),coef(ndim+1), &
    meth,1,omega,wksp(irpnt),iwksp(iipnt),iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfact2')
  return
!
!  propa = .false.
!
 20   call vicopy (maxnz,jcoef,iwksp(iipnt))
  maxt = maxnz - 1
  maxb = 0
  if (lvfill == 0) go to 26
  do 25 i = 1,lvfill
 25   call fills (maxt,iwksp(iipnt+1))
 26   nfactr = n*(maxt + 1)
  nfacti = maxt + 1
  call needw ('pfact2',1,iipnt,maxt**2,ier)
  if (ier < 0) return
  call needw ('pfact2',0,irpnt,nfactr+n,ier)
  if (ier < 0) return
!
  ifactr = irpnt
  ifacti = iipnt
  call vfill (nfactr,wksp(ifactr),0.0)
  do 40 j = 1,maxnz
     ip1 = ndim*(j - 1) + 1
     ip2 = n*(j - 1) + ifactr
     call vcopy (n,coef(ip1),wksp(ip2))
 40   continue
  irpnt = irpnt + nfactr
  iipnt = iipnt + maxt + 1
  call rowise (maxt+1,iwksp(ifacti),irwise)
  call needw ('pfact2',1,iipnt,maxt,ier)
  if (ier < 0) return
  if (ifctv == 0) then
    call icf(n,n,maxt,iwksp(ifacti+1),wksp(ifactr),wksp(ifactr+n), &
      meth,0,omega,wksp(irpnt),iwksp(iipnt),iflag)
  end if

  if (ifctv == 1) then
    call icfv(n,n,maxt,iwksp(ifacti+1),wksp(ifactr),wksp(ifactr+n), &
      meth,0,omega,wksp(irpnt),iwksp(iipnt),iflag)
  end if

  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfact2')
  return
end
subroutine pfact3 (coef,jcoef,wksp,iwksp,nn,meth,ier)
!
!*******************************************************************************
!
!! PFACT3 computes a point incomplete factorization.
!
!  Parameters:
!
!       n       order of system
!       meth    method of factorization
!                = 1   ic   (unmodified)
!                = 2   mic  (modified)
!       nfactr  amount of real workspace needed for factorization
!       ier     error flag
!
!  
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  integer jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
!
  n = nn
!
!  if requested, find out if matrix has property a.
!
  if (ipropa == 0) propa = .false.
  if (ipropa == 1) propa = .true.
  if (lvfill > 0) propa = .false.
  if (lvfill > 0) go to 20
  if (ipropa /= 2) go to 15
  call needw ('pfact3',1,iipnt,2*n,ier)
  if (ier < 0) return
  call prbndx (n,ndim,maxnz,jcoef,coef,iwksp(iipnt),iwksp(iipnt+n),propa,3)
  if (propa) ipropa = 1
  if (.not. propa) ipropa = 0
!
 15   if (.not. propa) go to 20
!
!  propa = .true.
!
  ifactr = irpnt
  nfactr = n
  nfacti = 0
  call needw ('pfact3',0,irpnt,nfactr+n,ier)
  if (ier < 0) return
  call rowise (maxnz,jcoef,irwise)
  call needw ('pfact3',1,iipnt,maxt*maxb,ier)
  if (ier < 0) return
  call vfill (n,wksp(ifactr),0.0)
  call vcopy (n,coef,wksp(ifactr))
  irpnt = irpnt + nfactr
  maxtp1 = maxt + 1
  call icfn (ndim,n,maxt,maxb,jcoef(2),jcoef(maxt+2), &
    wksp(ifactr),coef(ndim+1),coef(ndim*maxtp1+1), &
    meth,1,omega,wksp(irpnt),iwksp(iipnt),iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfact3')
  return
!
!  propa = .false.
!
 20   call vicopy (maxnz,jcoef,iwksp(iipnt))
  maxz = maxnz
  if (lvfill == 0) go to 26
  do 25 i = 1,lvfill
 25   call filln (maxz,iwksp(iipnt))
 26   nfactr = n*maxz
  nfacti = maxz
  call needw ('pfact3',1,iipnt,maxz,ier)
  if (ier < 0) return
  call needw ('pfact3',0,irpnt,nfactr,ier)
  if (ier < 0) return
!
  ifactr = irpnt
  ifacti = iipnt
  call vfill (nfactr,wksp(ifactr),0.0)
  do 40 j = 1,maxnz
     ip1 = ndim*(j - 1) + 1
     ip2 = n*(j - 1) + ifactr
     call vcopy (n,coef(ip1),wksp(ip2))
 40   continue
  irpnt = irpnt + nfactr
  iipnt = iipnt + maxz
  call rowise (maxz,iwksp(ifacti),irwise)
  call needw ('pfact3',0,irpnt,n,ier)
  if (ier < 0) return
  call move2 (n,n,maxz,iwksp(ifacti),wksp(ifactr), wksp(irpnt),iwksp(iipnt), &
    maxt,maxb)
  call needw ('pfact3',1,iipnt,maxt*maxb,ier)
  if (ier < 0) return
  ipt1 = ifacti + maxt + 1
  ipt2 = ifactr + n*(maxt + 1)
  call icfn (n,n,maxt,maxb,iwksp(ifacti+1),iwksp(ipt1),wksp(ifactr), &
    wksp(ifactr+n),wksp(ipt2),meth,0,omega,wksp(irpnt),iwksp(iipnt),iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfact3')
  return
end
subroutine pfactc (coef,jcoef,wksp,iwksp,nn,methh,ier)
!
!*******************************************************************************
!
!! PFACTC computes a point incomplete factorization.  (multicolor ordering)
!
!
!  Parameters:
!
!       n       order of system
!       meth    method of factorization
!                = 1   ic   (unmodified)
!                = 2   mic  (modified)
!       ier     error flag
!
!  
!
!
  integer jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
  n = nn
  meth = methh
!
!  if requested, find out if matrix has property a.
!
  if (ipropa == 0) propa = .false.
  if (ipropa == 1) propa = .true.
  if (ipropa /= 2) go to 15
  call needw ('pfactc',1,iipnt,2*n,ier)
  if (ier < 0) return
  call prbndx (n,ndim,maxnz,jcoef,coef,iwksp(iipnt),iwksp(iipnt+n),propa,1)
  if (propa) ipropa = 1
  if (.not. propa) ipropa = 0
!
 15   if (.not. propa) go to 30
!
!  propa = .true.
!
  ifactr = irpnt
  nfactr = n
  nfacti = 0
  call needw ('pfactc',0,irpnt,nfactr,ier)
  if (ier < 0) return
  call vcopy (n,coef,wksp(ifactr))
  irpnt = irpnt + nfactr
  ip1 = ndim + 1
  maxc = maxnz - 1
  call icfcp (ndim,ndim,n,maxc,jcoef(ip1),wksp(ifactr),coef(ip1),ncolor, &
    iwksp(ndt),iwksp(ndb),meth,1,iwksp(ipt),omega,iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfactc')
  return
!
!  propa = .false.
!
 30   ifactr = irpnt
  nfactr = n*maxnz
  nfacti = 0
  call needw ('pfactc',0,irpnt,nfactr,ier)
  if (ier < 0) return
  call vfill (nfactr,wksp(ifactr),0.0)
  do 45 j = 1,maxnz
     ip1 = ndim*(j - 1) + 1
     ip2 = n*(j - 1) + ifactr
     call vcopy (n,coef(ip1),wksp(ip2))
 45   continue
  irpnt = irpnt + nfactr
  ip1 = ndim + 1
  ip3 = ifactr + n
  maxc = maxnz - 1
  call icfcp (n,ndim,n,maxc,jcoef(ip1),wksp(ifactr),wksp(ip3),ncolor, &
    iwksp(ndt),iwksp(ndb),meth,0,iwksp(ipt),omega,iflag)
  if (iflag == 1) ier = -12
  if (iflag == 2) ier = 5
  if (iflag == 0) return
  call ershow (ier,'pfactc')
  return
end
subroutine pfs (n,b1,b2,x)
!
!*******************************************************************************
!
!! PFS does a penta-diagonal forward substitution.
!
!
!  This has the form (i+b1+b2)*x = y.
!
!  b1 and b2 are the first and second sub-diagonals.
!
!  Parameters:
!
!          n      order of system
!          b1     vector of length n-1 containing the first
!                  sub-diagonal elements
!          b2     vector of length n-2 containing the second
!                  sub-diagonal elements
!          x      on input, x contains y
!                 on output, x contains the solution to
!                  (i + b1 + b2)*x = y
!
!  
!
  dimension b1(1), b2(1), x(2)
!
  x(2) = x(2) - b1(1)*x(1)
  do 10 i = 3,n
 10   x(i) = x(i) - b1(i-1)*x(i-1) - b2(i-2)*x(i-2)
  return
end
subroutine pfsm (nn,nsize,b1,b2,x)
!
!*******************************************************************************
!
!! PFSM does a penta-diagonal forward substitution.
!
!
!  This has the form (i+b1+b2)*x = y
!     where b1 and b2 are subdiagonals of a system composed of
!     independent subsystems of size nsize.
!
!  Parameters:
!
!          n      order of system
!          nsize  order of the individual subsystems
!          b1     linear array of length n-1 containing the first
!                  sub-diagonal elements of the factorizations
!          b2     linear array of length n-2 containing the second
!                  sub-diagonal elements of the factorizations
!          x      on input, x contains y
!                 on output, x contains
!                 the solution to (i + b1 + b2)*x = y
!
!  
!
  dimension b1(nsize,1), b2(nsize,1), x(nsize,1)
!
  n = nn
  nsys = n/nsize
  do 10 j = 1,nsys
 10   x(2,j) = x(2,j) - b1(1,j)*x(1,j)
  do 20 i = 3,nsize
     do 15 j = 1,nsys
 15      x(i,j) = x(i,j) - b1(i-1,j)*x(i-1,j) - b2(i-2,j)*x(i-2,j)
 20   continue
  return
end
subroutine pgen (nn,p,ip,nc,ncolor)
!
!*******************************************************************************
!
!! PGEN constructs the permutation and its inverse for a given coloring.
!
!
!  Parameters:
!
!         n         order of system (= nn)
!         p         vector from prbndx upon input
!                   permutation vector upon output
!         ip        integer workspace vector upon input
!                   inverse permutation vector upon output
!         nc        number of points for each color (output)
!         ncolor    number of colors
!
!  
!
  integer p(1), ip(1), nc(1)
!
  n = nn
!
!  determine number of colors and number of elements for each
!     color.
!
  ncolor = 0
  do 5 i = 1,n
 5    nc(i) = 0
  do 10 i = 1,n
     ic = p(i)
     if (ncolor < ic) ncolor = ic
     nc(ic) = nc(ic) + 1
 10   continue
!
!  construct permutation vector.
!
  ip(1) = 1
  do 15 i = 2,ncolor
     ip(i) = ip(i-1) + nc(i-1)
 15   continue
  do 20 i = 1,n
     ic = p(i)
     p(i) = ip(ic)
     ip(ic) = ip(ic) + 1
 20   continue
!
!  construct inverse permutation vector.
!
  do 25 i = 1,n
     j = p(i)
     ip(j) = i
 25   continue
  return
end
subroutine pjac (diag,nn,r,z)
!
!*******************************************************************************
!
!! PJAC does the point Jacobi preconditioning.
!
!  Parameters:
!
!         diag    vector of length n containing the diagonal
!                  elements of the coefficient matrix
!         n       order of system (= nn)
!         r       residual
!         z       output vector
!
!  
!
  dimension r(1), z(1), diag(1)
!
  n = nn
  do 10 i = 1,n
 10   z(i) = r(i)/diag(i)
  return
end
subroutine pmdg (ndim,mdim,nn,maxnz,jcoef,coef,ncol,nc,p,ip,maxd,maxnew, &
  jcnew,wksp,iwksp,isym,ier)
!
!*******************************************************************************
!
!! PMDG permutes the matrix according to and index vector.
!
!
!     If room allows, it stores the permuted matrix in
!     diagonal format.  there will be enough room if the number
!     of diagonals needed does not exceed mdim.
!
!  Parameters:
!
!        ndim      row dimension of coef and jcoef arrays
!                   in defining routine
!        mdim      column dimension of coef and jcoef arrays in
!                   defining routine
!        n         order of system (active row size of coef and jcoef)
!        maxnz     active column size of coef and jcoef
!        jcoef     integer array of column numbers
!        coef      real array of coefficients
!        ncolor    number of colors in the permutation (= ncol)
!        nc        integer vector of length ncolor giving the
!                   number of nodes for each color
!        p         permutation vector
!        ip        inverse permuation vector
!        maxd      active columns in permuted matrix
!        maxnew    integer vector giving the number of diagonals
!                   created for each color
!        jcnew     integer array of size ncolor*max(maxnew(i))
!                   giving the diagonal numbers for each color
!        wksp      real workspace of length n
!        iwksp     integer workspace of length 2*n
!        isym      symmetric storage switch
!                   = 0    symmetric storage
!                   = 1    nonsymmetric storage
!        ier       error flag
!                  =  0   no errors detected
!                  = -9   mdim is less than the number of columns
!                          needed in coef to store the permuted
!                          matrix in diagonal format
!
!  
!
  integer jcoef(2), nc(1), p(1), maxnew(1), jcnew(ncol,1), iwksp(1), ip(1)
  dimension coef(ndim,1), wksp(1)
!
!
  n = nn
  ncolor = ncol
!
!  fill out rest of matrix if symmetric storage is used.
!
  if (isym /= 0) go to 2
  maxd = 2*maxnz - 1
  if (mdim < maxd) ier = -9
  if (ier < 0) return
!
  do 1 j = 2,maxnz
     ind = jcoef(j)
     len = n - ind
     jcol = maxnz + j - 1
     jcoef(jcol) = -ind
     call vfill (ind,coef(1,jcol),0.0)
     call vcopy (len,coef(1,j),coef(ind+1,jcol))
 1    continue
  maxnz = maxd
!
!  determine the number of created diagonals.
!
 2    do 5 i = 1,ncolor
     maxnew(i) = 1
     jcnew(i,1) = 0
 5    continue
  do 35 j = 2,maxnz
     ind = jcoef(j)
     do 10 i = 1,n
        iwksp(n+i) = i + ind
        if (coef(i,j) == 0.0) iwksp(n+i) = i
 10      continue
     call vscati (n,iwksp(n+1),p,iwksp)
     call vgathi (n,p,iwksp,iwksp)
     do 15 i = 1,n
 15      iwksp(i) = iwksp(i) - i
     ist = 1
     do 30 k = 1,ncolor
        ncc = nc(k)
        ied = ist + ncc - 1
        lim = maxnew(k)
        do 25 i = ist,ied
           id = iwksp(i)
           do 20 jj = 1,lim
              if (jcnew(k,jj) == id) go to 25
 20            continue
           lim = lim + 1
           maxnew(k) = lim
           if (lim > mdim) go to 40
           jcnew(k,lim) = id
 25         continue
        ist = ist + ncc
 30      continue
 35   continue
!
!  determine maxd.
!
 40   maxd = -1
  do 45 k = 1,ncolor
 45   maxd = max (maxd,maxnew(k))
  if (mdim < maxd) ier = -9
  if (ier < 0) return
!
!  permute matrix.
!
  do 55 j = 1,maxnz
     do 50 i = 1,n
 50      wksp(i) = coef(i,j)
     call vscatr (n,wksp,p,coef(1,j))
 55   continue
!
!  rearrange rows.
!
  ist = 1
  do 85 k = 1,ncolor
     ncc = nc(k)
     ied = ist + ncc - 1
     lim = maxnew(k)
     do 62 l = 1,lim
        jcol = jcnew(k,l)
        iwksp(n+jcol) = l
 62      continue
     do 80 i = ist,ied
        iip = ip(i)
        do 60 j = 2,maxnz
 60         wksp(j) = coef(i,j)
        do 63 j = 2,maxd
 63         coef(i,j) = 0.0
        do 75 j = 2,maxnz
           if (wksp(j) == 0.0) go to 75
           icol = p(iip + jcoef(j)) - i
           l = iwksp(n+icol)
           coef(i,l) = wksp(j)
 75         continue
 80      continue
     ist = ist + ncc
 85   continue
  return
end
subroutine pneu (suba,coef,jcoef,wfac,jwfac,d,nd,wksp,nn,r,z)
!
!*******************************************************************************
!
!! PNEU computes a point Neumann polynomial approximation to inverse(A).  
!
!  the output vector is --
!         z = p  (a)*r
!              np
!
!  Parameters:
!
!         suba    matrix-vector multiplication routine
!         d       vector of length n giving the diagonal elements
!                  of the matrix
!         nd      the degree of the polynomial desired
!         wksp    workspace of length n
!         n       order of system (= nn)
!         r       residual
!         z       output vector
!
!  
!
  external  suba
  dimension r(1), d(1), z(1), wksp(1), coef(1), jcoef(2),wfac(1), jwfac(1)
!
  n = nn
  do 10 i = 1,n
 10   z(i) = r(i)/d(i)
  if (nd <= 0) return
!
  do 20 k = 1,nd
     call suba (coef,jcoef,wfac,jwfac,n,z,wksp)
     do 15 i = 1,n
 15      z(i) = z(i) + (r(i) - wksp(i))/d(i)
 20   continue
  return
end
subroutine pointr (icall,wksp,iwksp,ier)
!
!*******************************************************************************
!
!! POINTR adjusts pointers according to IFACT.
!
!  Parameters:
!
!       icall     indicates beginning or ending call
!                  = 1  for beginning
!                  = 2  for ending
!       wksp      real workspace vector
!       iwksp     integer workspace vector
!
!  
!
  integer   iwksp(1)
  dimension wksp(1)
!
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
!
!  data common blocks
!
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  if (icall == 2) go to 15
!
!  initialize pointers.
!
  if (ifact == 0) return
  iipnt = 1
  irpnt = 1
  nfactr = 0
  nfacti = 0
  ifactr = 1
  ifacti = 1
  return
!
!  reset pointers for return
!
 15   if (ier < 0) return
  if (nfacti == 0) go to 20
  call vicopy (nfacti,iwksp(ifacti),iwksp)
  iipnt = nfacti + 1
  ifacti = 1
 20   if (nfactr == 0) return
  call vcopy (nfactr,wksp(ifactr),wksp)
  iwkpt2 = iwkpt2 - ifactr + 1
  irpnt = nfactr + 1
  ifactr = 1
  return
end
subroutine ppii (suba,coef,jcoef,wfac,jwfac,ainf,alpha,beta,nd,wksp,nn,r,z)
!
!*******************************************************************************
!
!! PPII computes the least squares polynomial approximation to inverse(A). 
!
!  the output vector is --
!         z = p  (a)*r
!              np
!
!  Parameters:
!
!         suba    matrix-vector multiplication routine
!         ainf    the infinity norm of matrix a
!         alpha,  the least squares weighting factors
!          beta
!         nd      the degree of the polynomial desired
!         wksp    workspace of length 2*n
!         n       order of system (= nn)
!         r       residual
!         z       output vector
!
!  
!
  external  suba
  dimension r(1), z(1), wksp(1), coef(1), jcoef(2),wfac(1), jwfac(1)
!
!
  n = nn
  np1 = n + 1
  al = alpha
  be = beta
!
  c1 = ((al+be+2.0)*(al+be+3.0))/(ainf*(al+2.0)*(al+be+2.0))
  do 10 i = 1,n
 10   z(i) = c1*r(i)
  if (nd <= 0) return
!
  do 15 i = 1,n
 15   wksp(i) = r(i)
  do 35 k = 1,nd
     fk = float (k)
     c1 = ((2.0*fk+al+be+2.0)*(2.0*fk+al+be+3.0))/ &
       (ainf*(fk+al+2.0)*(fk+al+be+2.0))
     c2 = (fk*(fk+be)*(2.0*fk+al+be))/ &
       ((fk+al+1.0)*(fk+al+be+1.0)*(2.0*fk+al+be+2.0))
     call suba (coef,jcoef,wfac,jwfac,n,z,wksp(np1))
     do 20 i = 1,n
 20      wksp(n+i) = r(i) - wksp(n+i)
     do 25 i = 1,n
 25      wksp(i) = wksp(i+n) + c2*wksp(i)
     do 30 i = 1,n
 30      z(i) = z(i) + c1*wksp(i)
 35   continue
  return
end
subroutine prbblk (ncol,ndis,iblock,lbhb,p,ip,propa)
!
!**********************************************************************
!
!! PRBBLK determines if the matrix has block property A.
!
!
!     see routine prbndx for an explanation of the algorithm
!     (block structure)
!
!  input parameters --
!
!         ncolor   number of diagonal blocks
!         ndis     number of distinct diagonal blocks
!         iblock   integer array of size 3 by ndis by max(lbhb(i))
!                   giving block constants
!         lbhb     integer vector of size ndis giving the number
!                   of diagonal blocks for each distinct block size.
!         p,ip     integer workspace vectors of length ncolor
!
!  output parameters --
!
!        p      contains information for constructing the permutation
!               array upon output
!        propa  a logical variable which is set to .true. if the
!               matrix has block property a and .false. otherwise
!
!
!****************************************************************************
!
!  
!
  integer   p(1), ip(1), iblock(3,ndis,1), lbhb(1)
  logical   propa
!
!  specifications for local variables
!
  integer   first, old, young, curtyp, type
!
!
     ncolor = ncol
     ndist = ndis
     index = 1
     do 5 i = 1,ncolor
        p(i) = 0
        ip(i) = 0
 5       continue
!
!  handle the first set of points until some adjacent points
!  are found
!
     first = 1
!
 10      p(first) = first
     if (ndist > 1) index = first
     maxnz = lbhb(index)
     if (maxnz > 1) go to 20
!
!  search for next entry that has not been marked
!
     if (first == ncolor) go to 65
     do 15 i = first+1,ncolor
        if (p(i) /= 0) go to 15
        first = i
        go to 10
 15      continue
     go to 65
!
!  first set of adjacent points found
!
 20      next = 1
     last = 1
     ip(1) = first
!
!  loop over labeled points indicated in the stack stored in
!  the array ip
!
 25      k = ip(next)
     curtyp = p(k)
     nxttyp = -curtyp
     if (ndist > 1) index = k
     maxnz = lbhb(index)
     if (maxnz <= 0) go to 55
     do 50 j = 1,maxnz
        jcol = k + iblock(1,index,j)
!
!  determine if element (k,j) is a diagonal element or zero.
!
        if (jcol < 1 .or. jcol > ncolor .or. jcol == k) go to 50
        if (iblock(3,index,j) == 0) go to 50
!
        type = p(jcol)
!
!     the following is a five way case statement dealing with the
!     labeling of the adjacent node.
!
!  case i.  if the adjacent node has already been labeled with
!              label equal to nxttyp, then skip to the next adjacent
!              node.
!
     if (type == nxttyp) go to 50
!
!  case ii.  if the adjacent node has not been labeled yet label
!               it with nxttyp and enter it in the stack
!
     if (type /= 0) go to 30
        last = last + 1
        ip(last) = jcol
        p(jcol) = nxttyp
        go to 50
!
!  case iii.  if the adjacent node has already been labeled with
!                opposite color and the same father seed, then there
!                is an irrecoverable color conflict.
!
 30      if (type == curtyp) go to 999
!
!  case iv.  if the adjacent node has the right color and a different
!               father node, then change all nodes of the youngest fathe
!               node to point to the oldest father seed and retain the
!               same colors.
!
     if (type * nxttyp < 1) go to  40
        old   = min ( iabs(type), iabs(nxttyp) )
        young = max ( iabs(type), iabs(nxttyp) )
        do  35 l = young,ncolor
           if (iabs(p(l)) == young) p(l) = isign(old, p(l))
 35         continue
        curtyp = p(k)
        nxttyp = -curtyp
        go to 50
!
!  case v.  if the adjacent node has the wrong color and a different
!              father node, then change all nodes of the youngest father
!              node to point to the oldest father node along with
!              changing their colors.  since until this time the
!              youngest father node tree has been independent no other
!              color conflicts will arise from this change.
!
 40      old   = min ( iabs(type), iabs(nxttyp) )
     young = max ( iabs(type), iabs(nxttyp) )
     do  45 l = young,ncolor
        if (iabs(p(l)) == young) p(l) = isign(old, -p(l))
 45      continue
     curtyp = p(k)
     nxttyp = -curtyp
!
!
!end of case statement
!
 50      continue
!
!  advance to next node in the stack
!
 55      next = next + 1
     if (next <= last) go to 25
!
!  all nodes in the stack have been removed
!
!  check for nodes not labeled.  if any are found
!  start the labeling process again at the first
!  node found that is not labeled.
!
     do 60 i = first+1,ncolor
        if (p(i) /= 0) go to 60
           first = i
           go to 10
 60      continue
!
!  all nodes are now typed either red or black.
!  red-black ordering possible.
!
 65      propa = .true.
     do 70 i = 1,ncolor
        if (p(i) >= 0) p(i) = 1
        if (p(i) <= 0) p(i) = 2
 70      continue
     return
!
!.... type conflict
!
 999  propa = .false.
  return
end
subroutine prbndx (nn,ndim,maxnzz,jcoef,coef,p,ip,propa,nstore)
!
!**********************************************************************
!
!! PRBNDX determines if the matrix has property A.
!
!
!     this algorithm assumes all neighbors of a particular node
!     are known.
!
!     (Purdue, diagonal data structures)
!     the algorithm is to mark the first node as red (arbitrary).
!     all of its adjacent nodes are marked black and placed in
!     a stack.  the remainder of the code pulls the first node
!     off the top of the stack and tries to type its adjacent nodes.
!     the typing of the adjacent point is a five way case statement
!     which is well commented below (see do loop 50).
!
!     the array p is used both to keep track of the color of a node
!     (red node is positive, black is negative) but also the father
!     node that caused the color marking of that point.  since
!     complete information on the adjacency structure is hard to come
!     by, this forms a link to enable the color change of a partial
!     tree when a recoverable color conflict occurs.
!
!     the array ip is used as a stack to point to the set of nodes
!     left to be typed that are known to be adjacent to the current
!     father node.
!
!
!*****************************************************************************
!
!  input parameters --
!
!        n      number of nodes.  (integer, scalar) (= nn)
!        ndim   row dimension of coef array
!        maxnz  maximum number of nonzeros per row
!        jcoef  integer data array
!        coef   real data array
!        p,ip   integer workspace vectors of length n
!        nstore data structure switch
!                = 1  Purdue
!                = 2  diagonal (symmetric or nonsymmetric)
!
!  output parameters --
!
!        p      contains information for constructing the permutation
!               array upon output
!        propa  a logical variable which is set to .true. if the
!               matrix has property a and .false. otherwise
!
  integer   p(1), ip(1), jcoef(ndim,1)
  dimension coef(ndim,1)
  logical   propa
!
  integer   first, old, young, curtyp, type
!
     n = nn
     maxnz = maxnzz
     do 5 i = 1,n
        p(i) = 0
        ip(i) = 0
 5       continue
!
!  handle the first set of points until some adjacent points
!  are found
!
     first = 1
!
 10      p(first) = first
     if (maxnz > 1) go to 20
!
!  search for next entry that has not been marked
!
     if (first == n) go to 65
     ibgn = first + 1
     do 15 i = ibgn,n
        if (p(i) /= 0) go to 15
        first = i
        go to 10
 15      continue
     go to 65
!
!  first set of adjacent points found
!
 20      next = 1
     last = 1
     ip(1) = first
!
!  loop over labeled points indicated in the stack stored in
!  the array ip
!
 25      k = ip(next)
     curtyp = p(k)
     nxttyp = -curtyp
     if (maxnz <= 0) go to 55
     do 50 j = 1,maxnz
        if (nstore == 1) jcol = jcoef(k,j)
        if (nstore >= 2) jcol = k + jcoef(j,1)
!
!  determine if element (k,j) is a diagonal element or zero.
!
        if (jcol < 1 .or. jcol > n .or. jcol == k) go to 50
        if (coef(k,j) == 0.0) go to 50
!
        type = p(jcol)
!
!     the following is a five way case statement dealing with the
!     labeling of the adjacent node.
!
!  case i.  if the adjacent node has already been labeled with
!              label equal to nxttyp, then skip to the next adjacent
!              node.
!
     if (type == nxttyp) go to 50
!
!  case ii.  if the adjacent node has not been labeled yet label
!               it with nxttyp and enter it in the stack
!
     if (type /= 0) go to 30
        last = last + 1
        ip(last) = jcol
        p(jcol) = nxttyp
        go to 50
!
!  case iii.  if the adjacent node has already been labeled with
!                opposite color and the same father seed, then there
!                is an irrecoverable color conflict.
!
 30      if (type == curtyp) go to 999
!
!  case iv.  if the adjacent node has the right color and a different
!               father node, then change all nodes of the youngest fathe
!               node to point to the oldest father seed and retain the
!               same colors.
!
     if (type * nxttyp < 1) go to  40
        old   = min ( iabs(type), iabs(nxttyp) )
        young = max ( iabs(type), iabs(nxttyp) )
        do  35 l = young,n
           if (iabs(p(l)) == young) p(l) = isign(old, p(l))
 35         continue
        curtyp = p(k)
        nxttyp = -curtyp
        go to 50
!
!  case v.  if the adjacent node has the wrong color and a different
!              father node, then change all nodes of the youngest father
!              node to point to the oldest father node along with
!              changing their colors.  since until this time the
!              youngest father node tree has been independent no other
!              color conflicts will arise from this change.
!
 40      old   = min ( iabs(type), iabs(nxttyp) )
     young = max ( iabs(type), iabs(nxttyp) )
     do  45 l = young,n
        if (iabs(p(l)) == young) p(l) = isign(old, -p(l))
 45      continue
     curtyp = p(k)
     nxttyp = -curtyp
!
!
!end of case statement
!
 50      continue
!
!  advance to next node in the stack
!
 55      next = next + 1
     if (next <= last) go to 25
!
!  all nodes in the stack have been removed
!
!  check for nodes not labeled.  if any are found
!  start the labeling process again at the first
!  node found that is not labeled.
!
     ibgn = first + 1
     do 60 i = ibgn,n
        if (p(i) /= 0) go to 60
           first = i
           go to 10
 60      continue
!
!  all nodes are now typed either red or black.
!  red-black ordering possible.
!
 65      propa = .true.
     do 70 i = 1,n
        if (p(i) >= 0) p(i) = 1
        if (p(i) <= 0) p(i) = 2
 70      continue
     return
!
!.... type conflict
!
 999  propa = .false.
  return
end
subroutine prep (coef,jcoef,wksp,iwksp,nn,nstore,ier)
!
!*******************************************************************************
!
!! PREP puts the diagonal entries of the matrix into column 1 of COEF.
!
!  Parameters:
!
!         n       dimension of matrix
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         wksp    workspace array of size n
!         iwksp   integer workspace
!         ier     error flag -- on return, values mean
!                      0 -- no errors detected
!                     -5 -- nonexistent diagonal element
!
!  
!
  integer   jcoef(2), iwksp(1)
  dimension coef(1), wksp(1)
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cmpart / mpstrt, mpart
!
  n = nn
  go to (5,10,10,15,15), nstore
 5    call prep1 (n,ndim,maxnz,jcoef,coef,ier)
  return
 10   call prep2 (n,ndim,maxnz,jcoef,coef,wksp,ier)
  return
 15   call needw ('prep',1,iipnt,2*n+1,ier)
  if (ier < 0) return
  call prep3 (n,maxnz,jcoef,jcoef(ndim+1),coef,mpart,iwksp,iwksp(n+2))
  mpstrt = iipnt
  iipnt = iipnt + mpart + 1
  return
end
subroutine prep1 (nn,ndim,maxnzz,jcoef,coef,ier)
!
!*******************************************************************************
!
!! PREP1 puts the diagonal elements of the matrix in column 1 of COEF  (Purdue data structure)
!
!  Parameters:
!
!         n       dimension of matrix ( = nn)
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array (= maxnzz)
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         ier     error flag -- on return, values mean
!                      0 -- no errors detected
!                     -5 -- nonexistent diagonal element
!
!  
!
  integer   jcoef(ndim,1)
  dimension coef(ndim,1)
!
  n = nn
  maxnz = maxnzz
!
  do 20 i = 1,n
     do 10 j = 1,maxnz
        if (jcoef(i,j) == i) go to 15
 10      continue
!
!  no diagonal entry for row i.
!
     ier = -5
     return
!
!  switch entries so that diagonal element is in column 1.
!
 15      if (j == 1) go to 20
     save = coef(i,j)
     coef(i,j) = coef(i,1)
     jcoef(i,j) = jcoef(i,1)
     coef(i,1) = save
     jcoef(i,1) = i
 20   continue
  return
end
subroutine prep2 (nn,ndim,maxnzz,jcoef,coef,wksp,ier)
!
!*******************************************************************************
!
!! PREP2 puts the diagonal entries of the matrix into column 1 of COEF.  (diagonal data structure)
!
!
!  Parameters:
!
!         n       dimension of matrix ( = nn)
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array (= maxnzz)
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         wksp    workspace array of size n
!         ier     error flag -- on return, values mean
!                      0 -- no errors detected
!                     -5 -- nonexistent diagonal element
!
!  
!
  integer   jcoef(2)
  dimension coef(ndim,1), wksp(1)
!
  n = nn
  maxnz = maxnzz
!
  do 10 j = 1,maxnz
     if (jcoef(j) == 0) go to 15
 10   continue
!
!  no main diagonal.
!
  ier = -5
  return
!
!  switch diagonals so that main diagonal is in column 1.
!
 15   if (j == 1) return
  do 20 i = 1,n
     wksp(i) = coef(i,1)
     coef(i,1) = coef(i,j)
     coef(i,j) = wksp(i)
 20   continue
  jcoef(j) = jcoef(1)
  jcoef(1) = 0
  return
end
subroutine prep3 (n,nz,ia,ja,a,m,np,iwksp)
!
!*******************************************************************************
!
!! PREP3 puts the diagonal elements of the matrix into the data structure.
!
!
!  It also adds duplicate triples, and defines the partition for matrix-vector
!     products.
!
!  Parameters:
!
!         n       number of equations
!         nz      length of ia, ja, and a vectors
!         ia      vector of i values
!         ja      vector of j values
!         a       vector of matrix coefficients
!         m       number of partitions (output)
!         np      on output, np contains the partition pointers.
!                  it must be at least m+1 in length.
!         iwksp   integer workspace vector of length n
!
!  
!
  integer ia(1), ja(1), iwksp(1), np(1)
  dimension a(1)
!
!  eliminate duplicates from the vectors by adding their
!     values in the a vector.  first, sort the vectors by
!     rows first and then by columns within each row.
!
  call vsrta1 (nz,ia,ja,a)
!
!  add duplicates.
!
  l = 1
  do 10 k = 2,nz
     i = ia(k)
     j = ja(k)
     aval = a(k)
     if (i == ia(l) .and. j == ja(l)) go to 5
     l = l + 1
     ia(l) = i
     ja(l) = j
     a(l) = aval
     go to 10
 5       a(l) = a(l) + aval
 10   continue
  nz = l
!
!  put main diagonal elements first.
!
  do 20 k = 1,nz
 15      i = ia(k)
     j = ja(k)
     if (i /= j) go to 20
     if (i == k) go to 20
     val = a(k)
     ia(k) = ia(i)
     ja(k) = ja(i)
     a(k) = a(i)
     ia(i) = i
     ja(i) = i
     a(i) = val
     go to 15
 20   continue
!
!  define partitions.
!
  kbgn = n + 1
  krep = kbgn
  mm = 1
  np(1) = 1
 25   mm = mm + 1
  np(mm) = kbgn
  do 30 i = 1,n
 30   iwksp(i) = 0
  nval = 0
  if (kbgn > nz) go to 50
  do 40 k = kbgn,nz
     i = ia(k)
     j = ja(k)
     if (iwksp(i) == 1 .or. iwksp(i) == 3 .or.iwksp(j) >= 2) go to 40
     nval = nval + 1
     iwksp(i) = iwksp(i) + 1
     iwksp(j) = iwksp(j) + 2
     if (k == krep) go to 35
     at = a(krep)
     it = ia(krep)
     jt = ja(krep)
     a(krep) = a(k)
     ia(krep) = i
     ja(krep) = j
     a(k) = at
     ia(k) = it
     ja(k) = jt
 35      krep = krep + 1
     if (nval >= n) go to 45
 40   continue
 45   kbgn = krep
  go to 25
 50   m = mm - 1
  return
end
subroutine prich (nn,r,z)
!
!*******************************************************************************
!
!! PRICH does the Richardson preconditioning.
!
!  Parameters:
!
!         n       order of system (= nn)
!         r       residual
!         z       output vector
!
!  
!
  dimension r(1), z(1)
!
  n = nn
  do 10 i = 1,n
 10   z(i) = r(i)
  return
end
subroutine psoln (nn,d,t1,t2,b1,b2,y,x)
!
!*******************************************************************************
!
!! PSOLN solves the system A*x = y for x, where A is a pentadiagonal system.  
!
!  d, t1, t2, b1, and b2 contain
!     the main, first and second super, and first and second sub
!     diagonals, respectively, of the factorization.
!
!  Parameters:
!
!          n      order of system
!          d      vector of length n containing the diagonal
!                  elements of the factorization matrix
!          t1     vector of length n-1 containing the first
!                  super-diagonal elements of the factorization
!          t2     vector of length n-2 containing the second
!                  super-diagonal elements of the factorization
!          b1     vector of length n-1 containing the first
!                  sub-diagonal elements of the factorization
!          b2     vector of length n-2 containing the second
!                  sub-diagonal elements of the factorization
!          y      the right-hand side
!          x      the solution to ax = y
!
!  
!
  dimension d(1), t1(1), t2(1), b1(1), b2(1), x(1), y(1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call pfs (n,b1,b2,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call pbs (n,t1,t2,x)
  return
end
subroutine psolnm (nn,nsize,d,t1,t2,b1,b2,y,x)
!
!*******************************************************************************
!
!! PSOLNM solves the system A*x = y for x, where a contains multiple pentadiagonal systems.  
!
!  d, t1, t2, b1, and b2 are
!     the main, first and second super, and the first and second
!     sub diagonals, respectively, of the factorization.
!
!  Parameters:
!
!          n      order of system
!          nsize  size of an individual subsystem
!          d      vector of length n containing the diagonal
!                  elements of the factorization matrix
!          t1     vector of length n-1 containing the first
!                  super-diagonal elements of the factorization
!          t2     vector of length n-2 containing the second
!                  super-diagonal elements of the factorization
!          b1     vector of length n-1 containing the first
!                  sub-diagonal elements of the factorization
!          b2     vector of length n-2 containing the second
!                  sub-diagonal elements of the factorization
!          y      the right-hand side
!          x      the solution to ax = y
!
!  
!
  dimension d(1), t1(1), t2(1), b1(1), b2(1), x(1), y(1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call pfsm (n,nsize,b1,b2,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call pbsm (n,nsize,t1,t2,x)
  return
end
subroutine pstop (ncall,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs,r, &
  z,zt,wk,nw,ier)
!
!*******************************************************************************
!
!! PSTOP is an interface to the PSTOPG routine using NULLPL and NULLPR.
!
  dimension zt(1), z(1), r(1), u(1), ubar(1), rhs(1), wk(1)
  dimension coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subql, subqr, nullpl, nullpr
!
  call pstopg (ncall,suba,subql,subqr,nullpl,nullpr,coef,jcoef,wfac,jwfac, &
    n,u,ubar,rhs,r,z,zt,wk,nw,ier)

  return
end
subroutine pstopg (ncall,suba,subql,subqr,precl,precr,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,r,z,zt,wk,nw,ier)
!
!*******************************************************************************
!
!! PSTOPG computes one of the stopping tests.
!
!
!  PSTOPG determine if the
!     iterative method has converged to a solution within the
!     error tolerance, zeta.  the stopping tests are --
!
!    (1)  (emax/emin) * sqrt ( (r ,zt)/(rhs,inv(q)*rhs) )
!    (2)  ( 1.0/emin) * sqrt ( (zt,zt)/(u,u) )
!    (3)  (emax/emin) * sqrt ( (zt,zt)/(inv(q)*rhs,inv(q)*rhs) )
!    (4)                sqrt ( (zt,zt)/(inv(q)*rhs,inv(q)*rhs) )
!    (5)                sqrt ( (r ,r )/(rhs,rhs) )
!    (6)                sqrt ( (u-ubar,u-ubar)/(ubar,ubar) )
!    (7)  (emax/emin) * sqrt ( (r,z)/(rhs,inv(ql)*rhs) )
!    (8)  ( 1.0/emin) * sqrt ( (z,z)/(u,u) )
!    (9)  (emax/emin) * sqrt ( (z,z)/(inv(ql)*rhs,inv(ql)*rhs) )
!   (10)                sqrt ( (z,z)/(inv(ql)*rhs,inv(ql)*rhs) )
!
!  here, emax and emin are estimates of the 2-norm of the iteration
!     matrix and its inverse.
!
! key parameters --
!
! ncall: = 0 for first call to pstop by accelerator
!        < 0 for recalc of bnorms, in the case that a new prec has
!            been calc'ed
!        > 0 for a routine call to check the stopping test
!
! iplr : = 0 the left and right preconditioning matrices are the
!            identity
!        = 1 the right prec is the identity
!        = 2 the left prec is the identity
!        = 3 neither the left nor the right prec matrix is the
!            identity
!
! r: the residual of the original system, if rhave = .true.
! z :  ql**(-1) r, if zhave = .true.
! zt : qr**(-1) z, if zthave = .true.
!
!
! this routine is admittedly quite overdesigned.  the idea was to have
! a general routine which would calculate the needed inner products
! with the absolute least amount of work, by determining which inner
! products already exist.
!
  dimension zt(1), z(1), r(1), u(1), ubar(1), rhs(1), wk(1)
  dimension coef(1), jcoef(2), wfac(1), jwfac(1)
  external suba, subql, subqr, precl, precr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
  logical init, nufact
  logical ipl, ipr
  logical risz, riszt, ziszt
  logical rhav, zhav, zthav, rcalc, zcalc, ztcalc
  logical udhv, rdhv, rzhv, rzthv, zdhv, zzthv, ztdhv
  logical udcal, rdcal, rzcal, rztcal, zdcal, zztcal, ztdcal
!
  dimension idlarr(10), idrarr(10), needbn(10)
  data idlarr /1,3,3,3,1,0,1,2,2,2/
  data idrarr /3,3,3,3,1,0,2,2,2,2/
  data needbn /1,0,1,1,1,0,1,0,1,1/
!
!
  nwusd = 0
  halt = .false.
  tiny = 500.0*srelpr
!
! get flags to tell us if there is any prec on the left or right.
  ipl = iplr == 1 .or. iplr == 3
  ipr = iplr == 2 .or. iplr == 3
! find equivalences between r, z, zt.
  risz = .not. ipl
  ziszt = .not. ipr
  riszt = risz .and. ziszt
! decode ntest.
  idl = idlarr(ntest)
  idr = idrarr(ntest)
  idot = 1 + (idl-1) + (idr-1)*3
!
  init = ncall == 0
  nufact = ncall < 0
  if (.not. (init .or. nufact)) go to 900
!
!  Initialization section.
!
  iv1 = 1
  iv2 = iv1 + n
!
!  compute bnorms, as necessary.
!
  if (needbn(ntest)==0) go to 750
  idle = idl
  if (idle==3 .and. ziszt) idle = 2
  if (idle==2 .and. risz)  idle = 1
  idre = idr
  if (idre==3 .and. ziszt) idre = 2
  if (idre==2 .and. risz)  idre = 1
  idp  = 0
  idlp = 0
  idrp = 0
  nwusd = 0
  if (nwusd > nw) go to 999
! calc ql(inv)*rhs, if necess.
  if (max(idle,idre)>1 .and. ipl) then
    nwusd = nwusd + n
    if (nwusd > nw) go to 999
    idp = idp + 1
    call precl (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,rhs,wk(1+n*(idp-1)))
    if (idle > 1) idlp = idlp + 1
    if (idre > 1) idrp = idrp + 1
  end if
! calc qr(inv)*ql(inv)*rhs, if necess.
  if (max(idle,idre)>2) then
    nwusd = nwusd + n
    if (nwusd > nw) go to 999
    idp = idp + 1
    if (idp == 1)call precr (coef,jcoef,wfac,jwfac,n,subql,suba,subqr, &
      rhs,wk(1+n*(idp-1)))
    if (idp == 2)call precr (coef,jcoef,wfac,jwfac,n,subql,suba,subqr, &
                       wk(1+n*(idp-2)),wk(1+n*(idp-1)))
    if (idle > 2) idlp = idlp + 1
    if (idre > 2) idrp = idrp + 1
  end if
! get needed dot.
  if (init .or. (idlp/=0 .or. idrp/=0)) then
    bnorm1 = seldot(n,1+idlp,rhs,wk(1),wk(1+n),wk(1+2*n), &
      1+idrp,rhs,wk(1),wk(1+n),wk(1+2*n))
    if (bnorm1 < 0e0) go to 998
    bnorm1 = amax1 (srelpr,sqrt(bnorm1))
  end if
  if (idlp==0 .or. idrp==0) bnorm = bnorm1
!
!  get ubar norm, as necessary.
!
 750  if (nufact) go to 900
  ubarnm = srelpr
  if (ntest == 6) ubarnm = sqrt(vdot (n,ubar,ubar))
!
!end of initialization phase.
!
!  now begin the actual stopping test section.
!
! notes on the strategy of this routine:
!     basically, what we're after in order to perform the stopping
! tests is certain dot products.  the needed dot products may already
! be available from the accelerator (in variables rrot, etc., as
! indicated in flags rdhav, etc.) otherwise, it will be necessary to
! compute these from the appropriate vectors.  these vectors in turn
! may already exist (in variables r, z, zt, as indicated b
! rhave, zhave, zthave), or it may be necessary to compute them.
! if they are computed by pstop, then the workspace is used to store
! them. furthermore, there are dependencies between the vectors: zt
! requires z, z requires r.  add to this the further complication that
! it may be possible to c optimize: if there is no left preconditioner,
! then r equals z, and so forth.
!     this routine attempts to get the necessary data to do the
! stopping test in the most optimal way.
!     a few notes on the semantics of variables.  the flag rhave tells
! whether the variable named r actually contains the residual; the flag
! rhav tells whether the residual exists somewhere - whether in r, z,
! zt or workspace.  if c nonzero, the variable ir tells where in the
! workspace the residual is (if it is in the workspace).  now, the
! variable rdhav tells whether rdot actually contains the dot product
! of r with itself.  unlike r and rhave, rdot and rdhav will actually
! be updated by pstop if they are calculated herein, or if rdot can be
! found from some other dot.
!     the variable rcalp indicates whether r was somewhere in workspace
! after pstop did its work.  the accelerator would like to know this,
! since it may want to circumvent letting pstop do a vector calculation
! if it can do it more efficiently.
!     for the initialization call (ncall=0), there is a dry run of
! the stopping test.  that is, the flags rhave, rcalp, rrhave, etc.
! are set to what they would be set in an actual call, but no actual
! vector calculations are done.  this is necessary so that the
! accelerator can plan ahead and take action to circumvent pstop doing
! lengthly calculations - e.g., calculating the residual using an a
! mult when the accelerator could do it simply by doing a saxpy.
!
!
 900  continue
!
! make temporaries for dot haves (modify the actual dot haves only
! if ncall>0)
  udhv = udhav
  rdhv = rdhav
  rzhv = rzhav
  rzthv = rzthav
  zdhv = zdhav
  zzthv = zzthav
  ztdhv = ztdhav
!
! evaluate vector haves.
!
  rhav  = rhave  .or. (zhave.and.risz)  .or. (zthave.and.riszt)
  zhav  = zhave  .or. (rhave.and.risz)  .or. (zthave.and.ziszt)
  zthav = zthave .or. (rhave.and.riszt) .or. (zhave .and.ziszt)
!
! take note that there are no vectors in the workspace.
!
  ir = 0
  iz = 0
  izt = 0
!
  iwfree = 1
!
!  Calculate R.
!
! find dot needs.
!
 102  assign 105 to lbldn
  go to 1100
! calculate whatever dots we can.
!
 105  assign 110 to lbldc
  go to 1300
! find vector needs.
!
 110  assign 115 to lblvn
  go to 1200
! get r.
!
 115  if (.not. rcalc) go to 120
  ir = iwfree
  iwfree = iwfree + n
  nwusd = iwfree-1
  if (init .or. nufact) go to 116
  if (nwusd > nw) go to 999
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(ir))
  call vexopy (n,wk(ir),rhs,wk(ir),2)
 116  rhav = .true.
! revise vector haves.
  if (.not. risz) go to 111
  iz = ir
  zhav = .true.
 111  if (.not. riszt) go to 120
  izt = ir
  zthav = .true.
!
!  Calculate z.
!
! calculate dots.
!
 120  assign 125 to lbldc
  go to 1300
! revise vector needs.
 125  assign 126 to lblvn
  go to 1200
! get z.
!
 126  if (.not. zcalc) go to 130
  iz = iwfree
  iwfree = iwfree + n
  nwusd = iwfree-1
  if (init .or. nufact) go to 127
  if (nwusd > nw) go to 999
  if (rhave) call precl (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,r,wk(iz))
  if (ir /= 0) call precl (coef,jcoef,wfac,jwfac,n,subql,suba, &
     subqr,wk(ir),wk(iz))
 127  zhav = .true.
! revise vector haves.
  if (.not. risz) go to 121
  ir = iz
  rhav = .true.
 121  if (.not. ziszt) go to 130
  izt = iz
  zthav = .true.
!
!  Calculate zt.
!
! calculate dots.
 130  assign 135 to lbldc
  go to 1300
!  revise vector needs ..
 135  assign 136 to lblvn
  go to 1200
! get zt.
 136  if (.not. ztcalc) go to 150
  izt = iwfree
  iwfree = iwfree + n
  nwusd = iwfree-1
  if (init .or. nufact) go to 137
  if (nwusd > nw) go to 999
  if (zhave)  call precr (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,z,wk(izt))
  if ((.not. zhave) .and. (rhave .and. risz)) then
    call precr (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,r,wk(izt))
  end if
  if (iz /= 0) then
    call precr (coef,jcoef,wfac,jwfac,n,subql,suba,subqr,wk(iz),wk(izt))
  end if

 137  zthav = .true.
! revise vector haves.
!
  if (.not. riszt) go to 131
  ir = izt
  rhav = .true.
 131  if (.not. ziszt) go to 150
  iz = izt
  zhav = .true.
!
!  Take care of details before going on to perform the stopping test.
!
! calculate whatever dots we can.
!
 150  assign 151 to lbldc
  go to 1300
! save vector calculation needs.
!
 151  rcalp  = ir  /= 0
  zcalp  = iz  /= 0
  ztcalp = izt /= 0
! head home, if ncall<=0.
!
  if (init .or. nufact) go to 950
!
! save dot have temporaries, if ncall>0.
!
  udhav  = udhv
  rdhav  = rdhv
  rzhav  = rzhv
  rzthav = rzthv
  zdhav  = zdhv
  zzthav = zzthv
  ztdhav = ztdhv
!
! get (u-ubar,u-ubar)
  if (ntest /= 6) go to 45
  uedot= 0.0
  do 40 i = 1,n
 40   uedot = uedot + (u(i) - ubar(i))**2
!
!  Stopping test computation section.
!
! at this point, all the needed dot products have been computed, and
! we are to actually perform the stopping test.
!
 45   go to (51,52,53,54,55,56,57,58,59,60), ntest
!
!  test 1
!
 51   if (rztdot < -srelpr) go to 998
  top = emax * sqrt (abs(rztdot))
  bottom = emin * bnorm1
  go to 80
!
!  test 2
!
 52   top = sqrt (abs(ztdot))
  bottom = emin * udnm
  go to 80
!
!  test 3
!
 53   top = emax * sqrt (abs(ztdot))
  bottom = emin * bnorm1
  go to 80
!
!  test 4
!
 54   top = sqrt (abs(ztdot))
  bottom = bnorm1
  go to 80
!
!  test 5
!
 55   top = sqrt (abs(rdot))
  bottom = bnorm1
  go to 80
!
!  test 6
!
 56   top = sqrt (abs(uedot))
  bottom = ubarnm
  go to 80
!
!  test 7
!
 57   if (rzdot < -srelpr) go to 998
  top = emax * sqrt (abs(rzdot))
  bottom = emin * bnorm1
  go to 80
!
!  test 8
!
 58   top = sqrt (abs(zdot))
  bottom = emin * udnm
  go to 80
!
!  test 9
!
 59   top = emax * sqrt (abs(zdot))
  bottom = emin * bnorm1
  go to 80
!
!  test 10
!
 60   top = sqrt (abs(zdot))
  bottom = bnorm1
  go to 80
!
 80   if (bottom < tiny) bottom = tiny
  stptst = top/bottom
  call ckconv (ier)
  if (ier < 0) go to 950
  halt = .false.
  if (top < bottom*zeta) halt = .true.
!
! done with the stopping test, head home.
  go to 950
!
!  Section to calculate dot-needs.
!
! here, we consider which dot products the stopping test needs, and
! see whether the needed dot products are currently nonexistent and
! thus must be calculated.
!
 1100 continue
!
! spread any dot information to other dots, as possible.
  if (risz) then
    if (rdhv) then
      rzdot = rdot
      rzhv = .true.
      zdot = rdot
      zdhv = .true.
    end if
    if (rzhv) then
      rdot = rzdot
      rdhv = .true.
      zdot = rzdot
      zdhv = .true.
    end if
    if (zdhv) then
      rzdot = zdot
      rzhv = .true.
      rdot = zdot
      rdhv = .true.
    end if
  end if
!
  if (ziszt) then
    if (zdhv) then
      zztdot = zdot
      zzthv = .true.
      ztdot = zdot
      ztdhv = .true.
    end if
    if (zzthv) then
      zdot = zztdot
      zdhv = .true.
      ztdot = zztdot
      ztdhv = .true.
    end if
    if (ztdhv) then
      zdot = ztdot
      zdhv = .true.
      zztdot = ztdot
      zzthv = .true.
    end if
  end if
!
  if (riszt) then
    if (rdhv) then
      rztdot = rdot
      rzthv = .true.
      ztdot = rdot
      ztdhv = .true.
    end if
    if (rzthv) then
      rdot = rztdot
      rdhv = .true.
      ztdot = rztdot
      ztdhv = .true.
    end if
    if (ztdhv) then
      rztdot = ztdot
      rzthv = .true.
      rdot = ztdot
      rdhv = .true.
    end if
  end if
!
! figure out which dots actually need to be calculated.
 1103 udcal = (needbn(ntest)==0 .and. ntest/=6) .and. .not.udhv
  rdcal  =  idot==1                 .and. .not.rdhv
  rzcal  = (idot==2 .or. idot==4) .and. .not.rzhv
  rztcal = (idot==3 .or. idot==7) .and. .not.rzthv
  zdcal  =  idot==5                 .and. .not.zdhv
  zztcal = (idot==6 .or. idot==8) .and. .not.zzthv
  ztdcal =  idot==9                 .and. .not.ztdhv
  go to lbldn
!
!  Section to calculate vector-needs.
!
! here, we see which vectors have to be calculated in order to
! satisfy the dot calculation needs.
!
 1200 continue
  ztcalc = (rztcal.or.zztcal.or.ztdcal) .and. .not.zthav
  zcalc  = (rzcal .or.zdcal .or.zztcal .or. ztcalc) .and. .not.zhav
  rcalc  = (rdcal .or.rzcal .or.rztcal .or. zcalc) .and. .not.rhav
  go to lblvn
!
!  Dot product calculation section.
!
! here, we calculate whatever dot products can be calculated from the
! currently existing vectors.
!
! first locate where the needed vectors are.
 1300 if (rhave) locr = 1
  if (zhave .and. risz) locr = 2
  if (zthave .and. riszt) locr = 3
  if (ir /= 0) locr = 4
!
  if (rhave .and. risz) locz = 1
  if (zhave) locz = 2
  if (zthave .and. ziszt) locz = 3
  if (iz /= 0) locz = 4
!
  if (rhave .and. riszt) loczt = 1
  if (zhave .and. ziszt) loczt = 2
  if (zthave) loczt = 3
  if (izt /= 0) loczt = 4
!
! now calculate whatever dot products we can.
!
!********** get udnm.
  if (.not. udcal) go to 1350
  if ((in > 5) .and. (mod(in,5) /= 0)) go to 1350
  uold = udnm
  if (init .or. nufact) go to 1349
  udnm = sqrt ( abs ( vdot (n,u,u) ) )
!     if ((in > 5) .and. (abs (udnm-uold) < udnm*zeta))
!    a           is3 = 1
  if (udnm < srelpr) udnm = 1.0
 1349 udhv = .true.
  assign 1350 to lbldn
  go to 1100
!
!********** get rdot.
 1350 if ( .not. (rdcal .and. rhav)) go to 1360
  if (init .or. nufact) go to 1359
  rdot = seldot (n,locr,r,z,zt,wk(ir),locr,r,z,zt,wk(ir))
 1359 rdhv = .true.
  assign 1360 to lbldn
  go to 1100
!
!********** get rzdot.
 1360 if (.not. (rzcal .and. rhav .and. zhav)) go to 1370
  if (init .or. nufact) go to 1369
  rzdot = seldot (n,locr,r,z,zt,wk(ir),locz,r,z,zt,wk(iz))
 1369 rzhv = .true.
  assign 1370 to lbldn
  go to 1100
!
!********** get rztdot.
 1370 if (.not. (rztcal .and. rhav .and. zthav)) go to 1380
  if (init .or. nufact) go to 1379
  rztdot = seldot (n,locr,r,z,zt,wk(ir),loczt,r,z,zt,wk(izt))
 1379 rzthv = .true.
  assign 1380 to lbldn
  go to 1100
!
!********** get zdot.
 1380 if (.not. (zdcal .and. zhav)) go to 1390
  if (init .or. nufact) go to 1389
  zdot = seldot (n,locz,r,z,zt,wk(iz),locz,r,z,zt,wk(iz))
 1389 zdhv = .true.
  assign 1390 to lbldn
  go to 1100
!
!********** get zztdot.
 1390 if (.not. (zztcal .and. zhav .and. zthav)) go to 1400
  if (init .or. nufact) go to 1399
  zztdot = seldot (n,locz,r,z,zt,wk(iz),loczt,r,z,zt,wk(izt))
 1399 zzthv = .true.
  assign 1400 to lbldn
  go to 1100
!
!********** get ztdot.
 1400 if (.not. (ztdcal .and. zthav)) go to 1410
  if (init .or. nufact) go to 1409
  ztdot = seldot (n,loczt,r,z,zt,wk(izt),loczt,r,z,zt,wk(izt))
 1409 ztdhv = .true.
  assign 1410 to lbldn
  go to 1100
!
 1410 continue
  go to lbldc

 950  nw = nwusd
  return
!
! splitting matrix is not positive definite
!
 998  ier = -7
  call ershow (ier,'pstop')
  go to 950
!
! insuff. real wksp
!
 999  ier = -2
  call ershow (ier,'pstop')
  go to 950
end
subroutine pstops (nn,r,z,u,ubar,ier)
!
!*******************************************************************************
!
!! PSTOPS performs a test to see if the iterative method has converged.
!
!     (cg and si routines)
!
!     the stopping tests are --
!
!    (1)  (emax/emin) * sqrt ( (r ,zt)/(rhs,inv(q)*rhs) )
!    (2)  ( 1.0/emin) * sqrt ( (zt,zt)/(u,u) )
!    (3)  (emax/emin) * sqrt ( (zt,zt)/(inv(q)*rhs,inv(q)*rhs) )
!    (4)                sqrt ( (zt,zt)/(inv(q)*rhs,inv(q)*rhs) )
!    (5)                sqrt ( (r ,r )/(rhs,rhs) )
!    (6)                sqrt ( (u-ubar,u-ubar)/(ubar,ubar) )
!    (7)  (emax/emin) * sqrt ( (r,z)/(rhs,inv(ql)*rhs) )
!    (8)  ( 1.0/emin) * sqrt ( (z,z)/(u,u) )
!    (9)  (emax/emin) * sqrt ( (z,z)/(inv(ql)*rhs,inv(ql)*rhs) )
!   (10)                sqrt ( (z,z)/(inv(ql)*rhs,inv(ql)*rhs) )
!
!
!  Parameters:
!
!         n       order of system
!         r       residual vector
!         z       pseudo-residual vector
!         u       solution estimate
!         ier     error flag
!                  =  0   no errors detected
!                  = -7   splitting matrix is not positive definite
!
!  
!
  dimension r(1), z(1), u(1), ubar(1)
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
  logical q1
  save    q1
  n = nn
  halt = .false.
  tiny = 500.0*srelpr
  nteste = ntest
  if (ntest > 6) nteste = nteste - 6
!
  go to (10,20,30,40,50,60), nteste
!
!  test 1
!
 10   if (rzdot >= 0.0) go to 15
  ier = -7
  call ershow (ier,'pstops')
  return
 15   emaxl = emax
  eminl = emin
  if (eminl < tiny) eminl = tiny
  tl = emaxl*sqrt (rzdot)
  tr = eminl*bnorm1
  stptst = tl/tr
  if (tl < tr*zeta) halt = .true.
  return
!
!  test 2
!
!  special procedure for zeroth iteration
!
 20   if (in >= 1) go to 25
  q1 = .false.
  udnm = 1.0
  stptst = sqrt (rzdot)
  if (stptst < tiny) halt = .true.
  return
!
!  in >= 1
!
!  test if udnm needs to be recomputed.
!
 25   if (q1) go to 28
  if ((in > 5)  .and.  (mod(in,5) /= 0)) go to 28
     uold = udnm
     udnm = sqrt ( vdot (n,u,u) )
     if (udnm < tiny) udnm = 1.0
     if ((in > 5)  .and. (abs (udnm-uold) < udnm*zeta)) q1 = .true.
!
!  compute stopping test.
!
 28   eminl = emin
  if (eminl < tiny) eminl = tiny
  tl = sqrt ( vdot (n,z,z) )
  tr = udnm*eminl
  stptst = tl/tr
  if (tl < tr*zeta) halt = .true.
  return
!
!  test 3.
!
 30   emaxl = emax
  eminl = emin
  if (eminl < tiny) eminl = tiny
  tl = emaxl*sqrt ( vdot (n,z,z) )
  tr = eminl*bnorm1
  stptst = tl/tr
  if (tl < tr*zeta) halt = .true.
  return
!
!  test 4.
!
 40   tl = sqrt ( vdot (n,z,z) )
  tr = bnorm1
  stptst = tl/tr
  if (tl < tr*zeta) halt = .true.
  return
!
!  test 5.
!
 50   tl = sqrt ( vdot (n,r,r) )
  tr = bnorm
  stptst = tl/tr
  if (tl < tr*zeta) halt = .true.
  return
!
!  test 6.
!
 60   sum = 0.0
  do 65 i = 1,n
 65   sum = sum + (u(i) - ubar(i))**2
  tl = sqrt (sum)
  tr = ubarnm
  stptst = tl/tr
  if (tl < tr*zeta) halt = .true.
  return
end
subroutine pvec (n,nv,iv,s,s1,idotw,it,il,ir,vect,dots,ndc,betas,gamma, &
  gamize,svec,wk,ier)
!
!*******************************************************************************
!
!! PVEC performs generalized Gram-Schmidt on a collection of vectors.
!
!
! it is used to update the table of direction vectors for
! generalized conjugate gradient methods per-iteration.
! note that this routine was intended to be rather general,
! including block conjugate gradient methods.
!
! params --
! n      - size of the vectors
! nv     - the size of the p-vector table.
!          ie., the table contains p(it-1), p(it-2),...,p(it-nv).
! iv     - number of p-vector-like objects we are dragging along.
!          eg., if iv=3, then we may be computing p, ap and q(inv)ap.
! s      - the block size for block conjugate gradient methods.
! s1     - indicates how many of the old p-vectors are to be used to
!          orthogonalize the new p-vector.
! idotw  - indicates the bandwidth of the matrix used to calculate
!          the betas.
!          generally equals s1, but if the h-matrix is symmetric
!          then = 1.
! it     - iteration number.  this routine calculates p(it).
! il,ir  - integers between 1 and iv.  indicate whether p, ap or
!          q(inv)ap
!          is to be used to calculate the inner product for
!          orthogonality.
! vect   - the p-vector table.
! dots   - workspace for the dot products.
! ndc    - the number of dot products that have already been
!          computed by formit.
! betas  - workspace for the betas.
! gamma  - an s by s matrix containing the coefficients from applying
!          gram schmidt to p(it).
! gamize - flag to indicate whether gram schmidt is to be applied
!          after p(it) is calculated.
! svec   - input packet of vectors to the p-vector calculation
!          process.
! wk     - workspace. must be of size s.
! ier    - error code
!
! array structure and indexing functions --
!
! vect(n,s,nv,iv)    jv
! svec(n,s,iv)       isv
! dots(s,s,idotw,s1) id
! betas(s,s,s1)      ib
! gamma(s,s)         -
!
  integer   s1, s, idotw
  dimension vect(1), svec(1)
  dimension dots(1), betas(1)
  dimension gamma(s,s)
  dimension wk(1)
  logical   gamize
  common / itcom4 / keygs, srelpr, keyzer
!
! define the necessary indexing functions.
!
  jv(i,j,k,l) = 1 + (i-1) + n*((j-1) + s*(mod(k,nv) + nv*(l-1)))
  isv(i,j,k) = 1 + (i-1) + n*((j-1) + s*(k-1))
  id(i,j,k,l) = 1 + (i-1) + s*((j-1) + s*((k-l) + idotw*mod(k,s1)))
  ib(i,j,k) = 1 + (i-1) + s*((j-1) + s*mod(k,s1))
!
  ier = 0
!
!  handle first iteration.
!
  if (it == 0  .or.  s1 <= 0) go to 1000
!
!  now handle general iteration.
!
!  first, calculate dot products (p(it-1),p(i)).
!
  ibgn = max (it-idotw,0)
  iend = it - 1 - ndc
  if (ibgn > iend) go to 10
  do 2 i = ibgn,iend
  do 2 j = 1,s
  do 2 k = 1,s
 2    dots(id(j,k,it-1,i)) = vdot (n,vect(jv(1,j,it-1,il)),vect(jv(1,k,i,ir)))
!
!  next, form all the new betas.
!
 10   ibgn = max (it-s1,0)
  iend = it - 1
  do 3 i = ibgn,iend
  do 34 l = 1,s
  do 35 k = 1,s
  wk(k) = -vdot (n,vect(jv(1,k,i,il)),svec(isv(1,l,ir)))
  jbgn = max (i-idotw+1,it-s1,0)
  jend = i - 1
  if (jend < jbgn) go to 35
  do 4 j = jbgn,jend
  do 4 m = 1,s
 4    wk(k) = wk(k) - dots(id(k,m,i,j))*betas(ib(m,l,j))
 35   continue
  call vcopy (s*s,dots(id(1,1,i,i)),gamma)
  call gauss (s,s,gamma,wk(1),betas(ib(1,l,i)),ier)
  if (ier /= 0) go to 999
 34   continue
 3    continue
!
!  now, get new p vectors.
!
  do 37 m = 1,iv
  do 37 i = ibgn,iend
  do 37 l = 1,s
  do 6 k = 1,s
 6    call vtriad (n,svec(isv(1,l,m)),svec(isv(1,l,m)), &
                 betas(ib(k,l,i)),vect(jv(1,k,i,m)),1)
 37   continue
!
!  copy new vectors into the table.
!
 1000 do 168 m = 1,iv
 168  call vcopy (n*s,svec(isv(1,1,m)),vect(jv(1,1,it,m)))
!
!  now calculate gamma and orthogonalize the new block of p-vectors
!
  call vfill (s*s,gamma,0.0)
  do 881 i = 1,s
 881  gamma(i,i) = 1.0
  if (.not. gamize) return
  do 879 i = 1,s
  if (i == 1) go to 882
  do 883 j = 1,i-1
 883  wk(j) = vdot (n,vect(jv(1,1,j,it)),vect(jv(1,1,i,it)))
  do 884 j = 1,i-1
  do 885 m = 1,iv
 885  call vtriad (n,vect(jv(1,i,it,m)),vect(jv(1,i,it,m)), &
      -wk(j),vect(jv(1,j,it,m)),1)
  do 886 k = j,i-1
 886  gamma(j,i) = gamma(j,i) - gamma(j,k)*wk(k)
 884  continue
 882  vnorm = sqrt(vdot(n,vect(jv(1,i,it,1)),vect(jv(1,i,it,1))))
  if (abs(vnorm) < srelpr**2) go to 999
  do 888 m = 1,iv
 888  call vtriad (n,vect(jv(1,i,it,m)),xxx,1.0/vnorm,vect(jv(1,i,it,m)),2)
  do 887 j = 1,i
 887  gamma(j,i) = gamma(j,i)/vnorm
 879  continue
  return
!
!  error return.
!
 999  ier = -100
  return
end
subroutine qrupd (ndim,nnz,nind,c,s,ucnbar,ucn,u,b,ier)
!
!*******************************************************************************
!
!! QRUPD updates the QR factorization of a banded upper Hessenberg matrix.
!
!
! parameters --
! ndim   - the current size of the Hessenberg matrix
! nnz    - the actual number of nonzeros in the band of the
!          Hessenberg matrix.  obviously, must be <= than nind.
! nind   - the bandwidth of the Hessenberg matrix, as stored
! c,s    - arrays which hold the cosines and sines of all the
!          rotations that have been performed so far
! u      - the new rightmost column of the Hessenberg matrix,
!          which is to be rotated
! b      - the element of the Hessenberg matrix to be zapped
!          by the new rotation
! ucnbar - the element of the Hessenberg matrix that b is to be
! ucn    - rotated into the new value of ucnbar, after the rotation
!
  dimension c(1), s(1), u(1)
!
! note -- due to the fortran implementation on the cyber 205, it is
! necessary to make ucnbar an array rather than a scalar.
!
  dimension ucnbar(1)
!
!  define the usual indexing functions.
!
  indv(i) = 1 + mod(i,nind)
  indu(i) = 1 + mod(i,nind+1)
!
!  indu is used to index u.
!
  if (ndim <= 1) return
!
!  apply all the old rotations to the column.
!
  jbgn = max(1,ndim-nnz+1)
  jend = ndim - 2
  if (jend < jbgn) go to 3
  do 2 j = jbgn,jend
  u1 = c(indv(j))*u(indu(j)) + s(indv(j))*u(indu(j+1))
  u2 =-s(indv(j))*u(indu(j)) + c(indv(j))*u(indu(j+1))
  u(indu(j)) = u1
  u(indu(j+1)) = u2
 2    continue
 3    continue
!
!  now proceed to form the new  2-by-2 rotation matrix.
!
  ucnb = ucnbar(1)
  denom = sqrt(ucnb*ucnb+b*b)
  if (abs(ucnb) >= 1.0e-40) denom = sign(denom,ucnb)
  if (abs(denom) < 1.0e-40) go to 999
  c(indv(ndim-1)) = ucnb/denom
  s(indv(ndim-1)) = b/denom
!
!  now apply the new rotation.
!
  u1 = c(indv(ndim-1))*u(indu(ndim-1))+s(indv(ndim-1))*u(indu(ndim))
  u2 =-s(indv(ndim-1))*u(indu(ndim-1))+c(indv(ndim-1))*u(indu(ndim))
  u(indu(ndim-1)) = u1
  u(indu(ndim)) = u2
  ucn = c(indv(ndim-1))*ucnb + s(indv(ndim-1))*b
  return
 999  ier = -14
  return
end
subroutine redblk (ndim,n,maxnz,coef,jcoef,p,ip,nstore,iwksp,ier)
!
!*******************************************************************************
!
!! REDBLK determines if the matrix has property A.
!
!
!  Parameters:
!
!        n        problem size
!        nstore   storage mode
!                  = 1  Purdue format
!                  = 2  symmetric diagonal format
!                  = 3  nonsymmetric diagonal format
!                  = 4  symmetric sparse format
!                  = 5  nonsymmetric sparse format
!        iwksp    integer workspace vector of length n
!        ier      error code
!                  =  0   no errors detected
!                  = -8   matrix does not have property a
!
!  common blocks
!
  integer jcoef(2), p(1), ip(1), iwksp(1)
  dimension coef(1)
  logical           propal
!
  go to (5,5,5,10,10), nstore
 5    call prbndx (n,ndim,maxnz,jcoef,coef,p,ip,propal,nstore)
  go to 15
 10   call bicol (n,maxnz,jcoef,jcoef(ndim+1),p,ip,iwksp,propal)
 15   if (propal) ier = 0
  if (.not. propal) ier = -8
  if (propal) return
  call ershow (ier,'redblk')
  return
end
subroutine rich1 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RICH1 drives the Richardson preconditioner.
!
  external accel, suba8, suba9, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom4 / keygs, srelpr, keyzer
!
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + n
  call split (accel,suba8,suba9,copy,copy,copy,copy,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  if (keygs == 1) irpnt = irpnt - n
  return
end
subroutine rich2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RICH2 drives the Richardson preconditioner.
!
  external accel, suba1, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  call split (accel,suba1,suba1,copy,copy,copy,copy,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine rich3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RICH3 drives the Richardson preconditioner.
!
  external accel, suba4, suba5, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  call split (accel,suba4,suba5,copy,copy,copy,copy,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine rich4 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RICH4 drives the Richardson preconditioner.
!
  external accel, suba12, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom4 / keygs, srelpr, keyzer
!
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba12,suba12,copy,copy,copy,copy,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine rich5 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RICH5 drives the Richardson preconditioner.
!
  external accel, suba13, suba14, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom4 / keygs, srelpr, keyzer
!
  iwkpt1 = irpnt
  if (keygs == 1) irpnt = irpnt + 2*n
  call split (accel,suba13,suba14,copy,copy,copy,copy,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  if (keygs == 1) irpnt = irpnt - 2*n
  return
end
subroutine rowise (maxnz,jcoef,irwise)
!
!*******************************************************************************
!
!! ROWISE determines whether a row-wise or diagonal-wise algorithm should be used.
!
!  This is for ic and ssor splittings with
!     diagonal storage.  this routine should be called after
!     final factorization is computed.
!
!  Parameters:
!
!          maxnz  number of number of diagonals stored
!          jcoef  vector of diagonal numbers for factorization
!                  array or matrix
!          irwise has a value upon output of
!                  0   if diagonal-wise algorithm should be used
!                  1   if row-wise algorithm should be used
!
!  
!
  integer   jcoef(2)
!
!  use a rowwise algorithm if  2 <= /jcoef(j)/ <= maxd
!     some j.
!
  maxd = 10
!
  irwise = 0
  do 15 j = 1,maxnz
     jcol = iabs(jcoef(j))
     if (jcol <= 1 .or. jcol > maxd) go to 15
     irwise = 1
     return
 15   continue
  return
end
subroutine rowsum (lda,n,maxnzz,a,x,isym)
!
!*******************************************************************************
!
!! ROWSUM computes the row sum of a matrix.
!
!
!  Parameters:
!
!        lda     leading dimension of array a
!        n       active size of array a
!        maxnz   number of columns in array a
!        a       array of size n by maxnz
!        x       vector of length n containing the row
!                 sum of a upon output
!        isym    symmetry switch
!                 = 0  matrix is a banded symmetric matrix
!                       with the diagonal in column one
!                 = 1  matrix is nonsymmetric
!
!  
!
  dimension a(lda,1), x(1)
!
  maxnz = maxnzz
  do 10 i = 1,n
 10   x(i) = 0.0
  do 20 j = 1,maxnz
     do 15 i = 1,n
 15      x(i) = x(i) + a(i,j)
 20   continue
  if (isym == 1 .or. maxnz <= 1) return
  do 30 j = 2,maxnz
     do 25 i = j,n
 25      x(i) = x(i) + a(i-j+1,j)
 30   continue
  return
end
subroutine rs6 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RS6 drives the reduced system method (Purdue storage with red-black coloring).
!
  external accel, suba10, suba11, subq1, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
!
!  compute red-black rhs.
!
  nr = iwksp(nc)
  nb = n - nr
  call needw ('rs6',0,irpnt,2*n,ier)
  if (ier < 0) return
  irhs = irpnt
  irpnt = irpnt + nr
  call vfill (2*n,wksp(irhs),0.0)
  call rsbegp (n,nr,ndim,maxnz,jcoef,coef,wksp(irhs),rhs,wksp(irpnt))
  iwkpt1 = irpnt
  irpnt = irpnt + n + nb
  call split (accel,suba10,suba11,subq1,subq1,subq1,subq1,copy,copy,noadp, &
    coef,jcoef,nr,u,ubar,wksp(irhs),wksp,iwksp,iparm,rparm,ier)
  call rsendp (n,nr,ndim,maxnz,jcoef,coef,u,rhs,wksp(iwkpt1))
  irpnt = irpnt - 2*n
  return
end
subroutine rs7 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RS7 drives the reduced system method (diagonal storage with red-black coloring).
!
  external accel, suba6, suba7, subq76, subq77, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
!
  n = nn
  t1 = timer (dummy)
  if (ifact == 1) call mfact (coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
!
!  compute red-black rhs.
!
  nr = iwksp(nc)
  nb = n - nr
  call needw ('rs7',0,irpnt,n,ier)
  if (ier < 0) return
  irhs = irpnt
  irpnt = irpnt + nr
  call rsbegd (n,n,nr,ndim,iwksp(maxnew),ndt,ndb,iwksp(jcnew), &
    coef,wksp(irhs),rhs,wksp(ifactr),wksp(irpnt))
  iwkpt1 = irpnt
  irpnt = irpnt + nb
  call split (accel,suba6,suba7,subq76,subq77,subq76,subq77,copy,copy,noadp, &
    coef,jcoef,nr,u,ubar,wksp(irhs),wksp,iwksp,iparm,rparm,ier)
  call rsendd (n,n,nr,ndim,iwksp(maxnew),ndt,ndb,iwksp(jcnew), &
    coef,u,rhs,wksp(ifactr))
  irpnt = irpnt - n
  return
end
subroutine rsad (nn,nsize,nrr,ndim,maxnew,ndtt,ndbb,jcnew,coef,c,b,dfac,wksp)
!
!*******************************************************************************
!
!! RSAD computes C = ( DR - T * inverse(DB) * B ) * B.
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     diagonal storage
!
!  Parameters:
!
!        n          order of system
!        nsize      size of an individual subsystem (if multiple
!                    systems)
!        nr         order of the red subsystem
!        ndim       row dimension of coef array
!        maxnew     number of columns in coef array
!        ndt        number of upper diagonals in diagonal block
!        ndb        number of lower diagonals in diagonal block
!        coef       real data structure
!        b          vector of length n containing bb behind br
!        c          vector of length nr containing cr
!        dfac       vector of length (1+nt+nb)*n to contain
!                    factorization of diagonal block upon output
!        wksp       workspace vector of length nb
!
!  
!
  integer   jcnew(2,1), maxnew(2)
  dimension coef(ndim,2), b(1), c(1), dfac(1), wksp(1)
!
  n = nn
  nr = nrr
  ndt = ndtt
  ndb = ndbb
  nrp1 = nr + 1
  nb = n - nr
  maxd = 1 + ndt + ndb
  maxz = maxnew(1) - maxd
  max2 = maxnew(2) - maxd
!
!  cr = dr*br.
!
  if (ndt+ndb > 0) go to 15
  do 10 i = 1,nr
 10   c(i) = coef(i,1)*b(i)
  go to 20
 15   call bmuln (ndim,nr,ndt,ndb,coef,coef(1,2),coef(1,ndt+2),b,c)
!
!  wksp = b*br
!
 20   if (maxz*max2 == 0) return
  do 25 i = 1,nb
 25   wksp(i) = 0.0
  call vaddd (ndim,2,nb,nr,max2,coef(nrp1,maxd+1),jcnew(2,maxd+1),wksp,b,-nr)
!
!  wksp = inv(db)*wksp
!
  if (ndt+ndb > 0) go to 35
  do 30 i = 1,nb
 30   wksp(i) = wksp(i)*dfac(i+nr)
  go to 40
 35   call bdsol (n,nb,nsize,ndt,ndb,dfac(nrp1),wksp,wksp,1)
!
!  cr = cr - t*wksp
!
 40   call vsubd (ndim,2,nr,nb,maxz,coef(1,maxd+1),jcnew(1,maxd+1), c,wksp,nr)
  return
end
subroutine rsap (ndimm,n,nr,maxnz,jcoef,coef,b,c,wksp)
!
!*******************************************************************************
!
!! RSAP computes  C = ( DR - T * inverse(DB) * B ) * B.
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     Purdue format
!
!  Parameters:
!
!        ndim       row dimension of coef,jcoef arrays
!        n          order of total system
!        nr         order of red subsystem
!        maxnz      number of columns in coef,jcoef arrays
!        jcoef      integer array of matrix column numbers
!        coef       real array of matrix coefficients
!        b,c        vectors of length nr
!        wksp       workspace array of length n + nb
!
!  
!
  integer   jcoef(ndimm,2)
  dimension coef(ndimm,2), b(1), c(1), wksp(1)
!
  ndim = ndimm
  do 10 i = 1,nr
 10   c(i) = coef(i,1)*b(i)
  if (maxnz <= 1) return
  np1 = n + 1
  nb = n - nr
  nrp1 = nr + 1
  maxm1 = maxnz - 1
  do 15 i = 1,n
 15   wksp(i) = 0.0
  call vaddp (ndim,ndim,nb,maxm1,coef(nrp1,2),jcoef(nrp1,2),wksp(nrp1), &
    b,wksp(np1))
  do 20 i = nrp1,n
 20   wksp(i) = wksp(i)/coef(i,1)
  call vsubp (ndim,ndim,nr,maxm1,coef(1,2),jcoef(1,2),c,wksp,wksp)
  return
end
subroutine rsatd (nn,nsize,nrr,ndim,maxnew,ndtt,ndbb,jcnew,coef,c,b,dfac,wksp)
!
!*******************************************************************************
!
!! RSATD computes  C = ((dr**t) - (b**t)*(db**(-t))*(t**t))*b.
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     diagonal storage
!
!  Parameters:
!
!        n          order of system
!        nsize      size of an individual subsystem (if multiple
!                    systems)
!        nr         order of the red subsystem
!        ndim       row dimension of coef array
!        maxnew     number of columns in coef array
!        ndt        number of upper diagonals in diagonal block
!        ndb        number of lower diagonals in diagonal block
!        coef       real data structure
!        b          vector of length n containing bb behind br
!        c          vector of length nr containing cr
!        dfac       vector of length (1+nt+nb)*n to contain
!                    factorization of diagonal block upon output
!        wksp       workspace vector of length nb
!
!  
!
  integer   jcnew(2,1), maxnew(2)
  dimension coef(ndim,2), b(1), c(1), dfac(1), wksp(1)
!
  n = nn
  nr = nrr
  ndt = ndtt
  ndb = ndbb
  nrp1 = nr + 1
  nb = n - nr
  maxd = 1 + ndt + ndb
  maxz = maxnew(1) - maxd
  max2 = maxnew(2) - maxd
!
!  cr = (dr**t)*br.
!
  if (ndt+ndb > 0) go to 15
  do 10 i = 1,nr
 10   c(i) = coef(i,1)*b(i)
  go to 20
 15   call bmulnt (ndim,nr,ndt,ndb,coef,coef(1,2),coef(1,ndt+2),b,c)
!
!  wksp = (t**t)*br
!
 20   if (maxz*max2 == 0) return
  do 25 i = 1,nb
 25   wksp(i) = 0.0
  call vadddt (ndim,2,nr,nb,maxz,coef(1,maxd+1),jcnew(1,maxd+1),wksp,b,nr)
!
!  wksp = (db**(-t))*wksp
!
  if (ndt+ndb > 0) go to 35
  do 30 i = 1,nb
 30   wksp(i) = wksp(i)*dfac(i+nr)
  go to 40
 35   call bdsolt (n,nb,nsize,ndt,ndb,dfac(nrp1),wksp,wksp)
!
!  cr = cr - (b**t)*wksp
!
 40   call vsubdt (ndim,2,nb,nr,max2,coef(nrp1,maxd+1),jcnew(2,maxd+1),c, &
    wksp,-nr)
  return
end
subroutine rsatp (ndimm,n,nr,maxnz,jcoef,coef,b,c,wksp)
!
!*******************************************************************************
!
!! RSATP computes  C = (dr - (b**t)*inv(db)*(t**t))*b.
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     Purdue format
!
!  Parameters:
!
!        ndim       row dimension of coef,jcoef arrays
!        n          order of total system
!        nr         order of red subsystem
!        maxnz      number of columns in coef,jcoef arrays
!        jcoef      integer array of matrix column numbers
!        coef       real array of matrix coefficients
!        b,c        vectors of length nr
!        wksp       workspace array of length n + nb
!
!  
!
  integer   jcoef(ndimm,2)
  dimension coef(ndimm,2), b(1), c(1), wksp(1)
!
  ndim = ndimm
  do 10 i = 1,nr
 10   c(i) = coef(i,1)*b(i)
  if (maxnz <= 1) return
  np1 = n + 1
  nb = n - nr
  nrp1 = nr + 1
  maxm1 = maxnz - 1
  do 15 i = 1,n
 15   wksp(i) = 0.0
  call vaddpt (ndim,ndim,nr,maxm1,coef(1,2),jcoef(1,2),wksp,b,wksp)
  do 20 i = nrp1,n
 20   wksp(i) = wksp(i)/coef(i,1)
  call vsubpt (ndim,ndim,nb,maxm1,coef(nrp1,2),jcoef(nrp1,2),c, &
    wksp(nrp1),wksp(np1))
  return
end
subroutine rsbegd (nn,nsize,nrr,ndim,maxnew,ndtt,ndbb,jcnew,coef,c,b,dfac,wksp)
!
!*******************************************************************************
!
!! RSBEGD computes  CR = br - t*inv(db)*bb.
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     diagonal storage
!
!  Parameters:
!
!        n          order of system
!        nsize      size of an individual subsystem (if multiple
!                    systems)
!        nr         order of the red subsystem
!        ndim       row dimension of coef array
!        maxnew     number of columns in coef array
!        ndt        number of upper diagonals in diagonal block
!        ndb        number of lower diagonals in diagonal block
!        coef       real data structure
!        b          vector of length n containing bb behind br
!        c          vector of length nr containing cr
!        dfac       vector of length (1+nt+nb)*n containing
!                    factorization of diagonal block upon input
!        wksp       workspace vector of length nb
!
!  
!
  integer   jcnew(2,1), maxnew(2)
  dimension coef(ndim,2), b(1), c(1), dfac(1), wksp(1)
!
  n = nn
  nr = nrr
  ndt = ndtt
  ndb = ndbb
  nrp1 = nr + 1
  nb = n - nr
  maxd = 1 + ndt + ndb
!
!  compute cr.
!
  do 10 i = 1,nr
 10   c(i) = b(i)
  call bdsol (n,nb,nsize,ndt,ndb,dfac(nrp1),b(nrp1),wksp,1)
  maxm1 = maxnew(1) - maxd
  call vsubd (ndim,2,nr,nb,maxm1,coef(1,maxd+1),jcnew(1,maxd+1),c,wksp,nr)
  return
end
subroutine rsbegp (n,nr,ndim,maxnz,jcoef,coef,c,b,wksp)
!
!*******************************************************************************
!
!! RSBEGP computes  cr = br - t*inv(db)*bb.
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     Purdue storage
!
!  Parameters:
!
!        n          order of system
!        nr         order of the red subsystem
!        ndim       row dimension of coef array
!        maxnz      number of columns in coef array
!        jcoef      integer data structure
!        coef       real data structure
!        b          vector of length n containing bb behind br
!        c          vector of length nr containing cr
!        wksp       workspace vector of length n
!
!  
!
  integer   jcoef(ndim,2)
  dimension coef(ndim,2), b(1), c(1), wksp(1)
!
  nrp1 = nr + 1
  do 10 i = 1,nr
 10   c(i) = b(i)
  if (maxnz <= 1) return
  do 15 i = nrp1,n
 15   wksp(i) = b(i)/coef(i,1)
  maxm1 = maxnz - 1
  call vsubp (ndim,ndim,nr,maxm1,coef(1,2),jcoef(1,2),c,wksp,wksp)
  return
end
subroutine rsendd (nn,nsize,nrr,ndim,maxnew,ndtt,ndbb,jcnew,coef,x,b,dfac)
!
!*******************************************************************************
!
!! RSENDD computes  xb = inv(db)*(bb - b*xr).
!
!                 a  = ( dr   t )
!                       ( b   db )
!
!     diagonal storage
!
!  Parameters:
!
!        n          order of system
!        nsize      size of an individual subsystem (if multiple
!                    systems)
!        nr         order of the red subsystem
!        ndim       row dimension of coef array
!        maxnew     number of columns in coef array
!        ndt        number of upper diagonals in diagonal block
!        ndb        number of lower diagonals in diagonal block
!        coef       real data structure
!        x          vector of length n containing  xr, xb
!        b          vector of length n containing bb in the last
!                    nb locations
!        dfac       vector of length (1+nt+nb)*n containing
!                    factorization of diagonal block upon input
!
!  
!
  integer   jcnew(2,1), maxnew(2)
  dimension coef(ndim,2), x(1), b(1), dfac(1)
!
  n = nn
  nr = nrr
  ndt = ndtt
  ndb = ndbb
  nrp1 = nr + 1
  nb = n - nr
  maxd = 1 + ndt + ndb
!
!  compute xb.
!
  do 10 i = nrp1,n
 10   x(i) = b(i)
  max2 = maxnew(2) - maxd
  call vsubd (ndim,2,nb,nr,max2,coef(nrp1,maxd+1),jcnew(2,maxd+1),x(nrp1),x,-nr)
  call bdsol (n,nb,nsize,ndt,ndb,dfac(nrp1),x(nrp1),x(nrp1),1)
  return
end
subroutine rsendp (n,nr,ndim,maxnz,jcoef,coef,x,b,wksp)
!
!*******************************************************************************
!
!! RSENDP computes  xb = inv(db)*(bb - b*xr).
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     Purdue format
!
!  Parameters:
!
!        n          order of matrix
!        nr         order of red subsystem
!        ndim       row dimension of ah and jah arrays
!        maxnz      number of columns in coef and jcoef arrays
!        jcoef      integer array of column numbers
!        coef       real array of matrix coefficients
!        x          vector of length n containing  xr, xb
!        b          vector of length n containing bb in the last
!                    nb locations
!        wksp       workspace array of length nb
!
!  
!
  integer   jcoef(ndim,2)
  dimension coef(ndim,2), x(1), b(1), wksp(1)
!
  nrp1 = nr + 1
  nb = n - nr
  do 10 i = nrp1,n
 10   x(i) = b(i)
  if (maxnz <= 1) go to 15
  maxm1 = maxnz - 1
  call vsubp (ndim,ndim,nb,maxm1,coef(nrp1,2),jcoef(nrp1,2),x(nrp1),x,wksp)
 15   do 20 i = nrp1,n
 20   x(i) = x(i)/coef(i,1)
  return
end
subroutine rsmatd (ndim,nrr,nb,maxnew,jcnew,dr,ah,ak,db,maxrss,jcrs,rs, &
  maxlim,isym,ier)
!
!*******************************************************************************
!
!! RSMATD computes RS = dr - ah*inv(db)*ak.
!
!  A has been permuted to red-black form --
!
!                   * dr  ah *
!             a =   *        *
!                   * ak  db *
!
!     (diagonal storage)
!
!      dr is nr x nr        ah is nr x nb
!      ak is nb x nr        db is nb x nb
!
!  definition of parameters --
!
!         ndim          row dimension of ah and ak arrays
!         nr            number of red points
!         nb            number of black points
!         maxnew        integer vector of length 2 indicating number
!                        of diagonals stored in ah and ak,
!                        respectively.
!         jcnew         integer array of diagonal numbers
!         dr            vector of length nr
!         ah            array of size nr by (maxnew(1)-1)
!         ak            array of size nb by (maxnew(2)-1)
!         db            vector of length nb
!         maxrs         number of columns needed to store reduced
!                        system (output)
!         jcrs          diagonal numbers for rs (output)
!         rs            array to contain reduced system
!         maxlim        maximum column width to be allowed for rs
!         isym          symmetry switch for rs matrix
!                        = 0   store only upper half of rs
!                        = 1   store all of rs
!         ier           error code
!                        =  0     no errors detected
!                        = -2     maxlim < maxrs
!
!  
!
  integer maxnew(2), jcnew(2,1), jcrs(1)
  dimension db(1), ak(ndim,1), ah(ndim,1), dr(1), rs(nrr,1)
!
  nr = nrr
  maxrs = 1
  jcrs(1) = 0
  do 5 i = 1,nr
 5    rs(i,1) = dr(i)
  maxh = maxnew(1) - 1
  maxk = maxnew(2) - 1
  do 35 lh = 1,maxh
     i = jcnew(1,lh+1) - nr
     ia1 = max (1,1-i)
     ib1 = min (nr,nb-i)
     do 30 lk = 1,maxk
        k = jcnew(2,lk+1) + nr
        l = i + k
        if (l < 0 .and. isym == 0) go to 30
        do 10 ld = 1,maxrs
           if (jcrs(ld) == l) go to 20
 10         continue
        if (maxrs == maxlim) go to 999
        maxrs = maxrs + 1
        ld = maxrs
        jcrs(maxrs) = l
        do 15 ii = 1,nr
 15         rs(ii,maxrs) = 0.0
 20         ist = max (ia1,1-l)
        ied = min (ib1,nr-l)
        do 25 m = ist,ied
 25         rs(m,ld) = rs(m,ld) - ah(m,lh)*ak(m+i,lk)/db(m+i)
 30      continue
 35   continue
  maxrss = maxrs
  return
!
!  error exit -- maxlim too small.
!
 999  ier = -2
  return
end
subroutine rsmatp (ndim,nrr,maxnzz,jcoef,coef,maxrss,jcrs,rs,maxlim,wksp, &
  iwksp,ier)
!
!*******************************************************************************
!
!! RSMATP computes RS = dr - ah*inv(db)*ak.
!
!  A has been permuted to red-black form --
!
!                   * dr  ah *
!             a =   *        *
!                   * ak  db *
!
!     (Purdue storage)
!
!      dr is nr x nr        ah is nr x nb
!      ak is nb x nr        db is nb x nb
!
!  definition of parameters --
!
!         ndim          row dimension of coef and jcoef arrays
!         nr            number of red points
!         maxnz         number of columns in coef and jcoef
!         jcoef         array of column indices
!         coef          array of matrix coefficients
!         maxrs         number of columns needed to store reduced
!                        system (output)
!         jcrs          column numbers for rs (output)
!         rs            array to contain reduced system
!         maxlim        maximum column width to be allowed for rs
!         wksp          workspace of length 2*nr
!         iwksp         integer workspace of length nr
!         ier           error code
!                        =  0     no errors detected
!                        = -2     maxlim < maxrs
!
!  
!
  integer jcoef(ndim,1), jcrs(nrr,1), iwksp(1)
  dimension coef(ndim,1), rs(nrr,1), wksp(1)
!
  nr = nrr
  maxnz = maxnzz
  maxrs = 1
  do 5 i = 1,nr
     rs(i,1) = coef(i,1)
     jcrs(i,1) = i
 5    continue
  do 50 j = 2,maxnz
     call vgathr (nr,coef,jcoef(1,j),wksp)
     do 10 i = 1,nr
 10      wksp(i) = coef(i,j)/wksp(i)
     do 45 jj = 2,maxnz
        call vgathr (nr,coef(1,jj),jcoef(1,j),wksp(nr+1))
        call vgathi (nr,jcoef(1,jj),jcoef(1,j),iwksp)
        do 15 i = 1,nr
 15         wksp(nr+i) = wksp(i)*wksp(nr+i)
        do 40 i = 1,nr
           jcol = iwksp(i)
           term = wksp(nr+i)
           if (jcol > nr) go to 40
           do 20 jjj = 1,maxrs
              if (jcrs(i,jjj) /= jcol) go to 20
              rs(i,jjj) = rs(i,jjj) - term
              go to 40
 20            continue
           if (maxrs == 1) go to 30
           do 25 jjj = 2,maxrs
              if (jcrs(i,jjj) /= i) go to 25
              rs(i,jjj) = rs(i,jjj) - term
              jcrs(i,jjj) = jcol
              go to 40
 25            continue
 30            if (maxrs == maxlim) go to 999
           maxrs = maxrs + 1
           do 35 ii = 1,nr
              jcrs(ii,maxrs) = ii
              rs(ii,maxrs) = 0.0
 35            continue
           rs(i,maxrs) = -term
           jcrs(i,maxrs) = jcol
 40         continue
 45      continue
 50   continue
  maxrss = maxrs
  return
!
!  error exit -- maxlim too small.
!
 999  ier = -2
  return
end
subroutine rsnsp (precon,accel,ndimm,mdimm,nn,maxnzz,coef,jcoef,p,ip,u, &
  ubar,rhs,wksp,iwksp,nw,inw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! RSNSP is the driver for explicitly computed reduced systems.
!
!
!  Parameters:
!
!       precon    preconditioning module
!       accel     acceleration module
!       coef      real matrix data array
!       jcoef     integer matrix data array
!       n         input integer.  order of the system (= nn)
!       u         input/output vector.  on input, u contains the
!                 initial guess to the solution.  on output, it
!                 contains the latest estimate to the solution.
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       wksp      real workspace of length nw
!       iwksp     integer workspace of length inw
!       nw        length of wksp upon input, amount used upon
!                  output
!       inw       length of iwksp upon input, amount used upon
!                  output
!       iparm     integer vector of length 30.  allows user to
!                 specify some integer parameters which affect
!                 the method.
!       rparm     real vector of length 30.  allows user to
!                 specify some real parameters which affect
!                 the method.
!       ier       output integer.  error flag.
!
!  
!
  external  accel, precon
  integer   iparm(30), jcoef(2), p(1), ip(1), iwksp(1)
  dimension coef(1), rhs(1), u(1), ubar(1), rparm(30), wksp(1)
!
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  ier = 0
  ndim = ndimm
  mdim = mdimm
  n = nn
  maxnz = maxnzz
  lenr = nw
  leni = inw
  irmax = 0
  iimax = 0
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,1,ier)
  timfac = 0.0
  call pointr (1,wksp,iwksp,ier)
!
!  call preparatory routines.
!
!  remove zeros from jcoef for Purdue data structure.
!
  if (nstore == 1) call adjust (n,ndim,maxnz,jcoef,1)
  call prep (coef,jcoef,wksp(irpnt),iwksp(iipnt),n,nstore,ier)
  if (ier < 0) then
     call ershow (ier,'rsnsp')
     go to 20
  end if
!
!  eliminate penalty-method dirichlet points, if requested.
!
  ielim = iparm(24)
  tol = rparm(15)
  if (ielim == 1) call elim (n,jcoef,coef,rhs,wksp,iwksp,tol)
!
!  determine symmetry of matrix.
!
  if (nstore == 1 .and. isymm == 2) call detsym(ndim,maxnz,coef,jcoef,n,isymm)
!
!  form reduced system matrix.
!
  call rsprep (coef,jcoef,wksp,iwksp,n,rhs,u,ubar,p,ip,nr,irs,ijcrs, &
    irsrhs,ier)
!
!  scale matrix.
!
  call scale (wksp(irs),iwksp(ijcrs),wksp,1,nr,u,ubar,wksp(irsrhs),ier)
  if (ier < 0) go to 20
!
!  call iterative routine.
!
  call precon (accel,wksp(irs),iwksp(ijcrs),nr,u,ubar,wksp(irsrhs),wksp, &
    iwksp,iparm,rparm,ier)
!
!  unscale matrix.
!
  call scale (wksp(irs),iwksp(ijcrs),wksp,2,nr,u,ubar,wksp(irsrhs),ier)
!
!  restore to original system.
!
  call rspost (coef,jcoef,wksp,iwksp,n,rhs,u,ubar,p,ip,nr,irs,ijcrs,ier)
!
!  restore zeros to jcoef for Purdue data structure.
!
 20   if (nstore == 1) call adjust (n,ndim,maxnz,jcoef,2)
  t2 = timer (dummy)
  timtot = t2 - t1
  iparm(18) = ipropa
  iparm(23) = isymm
  rparm(13) = timfac
  rparm(14) = timtot
  call echall (n,iparm,rparm,2,1,ier)
!
  call pointr (2,wksp,iwksp,ier)
  nw = irmax
  inw = iimax
  maxnzz = maxnz
  return
end
subroutine rspost (coef,jcoef,wksp,iwksp,nn,rhs,u,ubar,p,ip,nrr,irs,ijcrs,ier)
!
!*******************************************************************************
!
!! RSPOST is the postprocessor for explicitly-computed reduced systems.
!
!  Parameters:
!
!       coef      real matrix data array
!       jcoef     integer matrix data array
!       n         input integer.  order of the system (= nn)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       u         current solution estimate
!       ubar      exact solution vector (if known)
!       nr        order of the reduced system upon input
!       irs       pointer into wksp for reduced system matrix
!       ijcrs     pointer into wksp for reduced system integer
!                  array
!       ier       output integer.  error flag.
!
!  
!
  integer jcoef(2), iwksp(1), p(1), ip(1)
  dimension coef(1), rhs(1), u(1), ubar(1), wksp(1)
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / rscons / ndimrs, mdimrs, maxzrs
!
  n = nn
  nr = nrr
  nb = n - nr
!
!  update constants.
!
  ndim = ndimrs
  mdim = mdimrs
  maxnz = maxzrs
  irpnt = irs
  iipnt = ijcrs
!
!  compute xb.
!
  call needw ('rspost',0,irpnt,nb,ier)
  if (ier < 0) return
  if (nstore == 1) call rsendp (n,nr,ndim,maxnz,jcoef,coef,u,rhs,wksp(irpnt))
  if (nstore >= 2) call rsxbd (n,nr,ndim,iwksp(maxnew),iwksp(jcnew),coef,u,rhs)
!
!  unpermute matrix.
!
  call permut (coef,jcoef,p,ip,wksp,iwksp,2,n,u,ubar,rhs,ier)
  if (ier < 0) return
  return
end
subroutine rsprep (coef,jcoef,wksp,iwksp,nn,rhs,u,ubar,p,ip,nrr,irs, &
  ijcrs,irsrhs,ier)
!
!*******************************************************************************
!
!! RSPREP is the preprocessor for explicitly-computed reduced systems.
!
!  Parameters:
!
!       coef      real matrix data array
!       jcoef     integer matrix data array
!       n         input integer.  order of the system (= nn)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       u         current solution estimate
!       ubar      exact solution vector (if known)
!       nr        order of the reduced system upon output
!       irs       pointer into wksp for reduced system matrix
!       ijcrs     pointer into wksp for reduced system integer
!                  array
!       irsrhs    pointer into wksp for reduced system rhs
!       ier       output integer.  error flag.
!
!  
!
  integer jcoef(2), iwksp(1), p(1), ip(1)
  dimension coef(1), rhs(1), u(1), ubar(1), wksp(1)
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / rscons / ndimrs, mdimrs, maxzrs
!
  n = nn
!
!  permute matrix.
!
  call permut (coef,jcoef,p,ip,wksp,iwksp,1,n,u,ubar,rhs,ier)
  if (ier < 0) return
!
!  form reduced system matrix.
!
  nr = iwksp(nc)
  nb = iwksp(nc+1)
  irs = irpnt
  ijcrs = iipnt
  length = lenr - irpnt + 1
  call vfill (length,wksp(irpnt),0.0)
  if (nstore >= 2) go to 30
!
!  Purdue storage.
!
  call needw ('rsprep',0,irpnt,3*nr,ier)
  if (ier < 0) return
  call needw ('rsprep',1,iipnt,2*nr,ier)
  if (ier < 0) return
  lim1 = (lenr - 2*nr - irpnt + 1)/nr
  lim2 = (leni - nr - iipnt + 1)/nr
  maxlim = min(lim1,lim2)
  ip1 = irpnt + nr*maxlim
  ip2 = iipnt + nr*maxlim
  call rsmatp (ndim,nr,maxnz,jcoef,coef,maxrs,iwksp(ijcrs), &
    wksp(irs),maxlim,wksp(ip1),iwksp(ip2),ier)
  if (ier < 0) then
     call ershow (ier,'rsprep')
     return
  end if
  irpnt = irpnt + nr*maxrs
  iipnt = iipnt + nr*maxrs
  go to 45
!
!  diagonal storage.
!
 30   call needw ('rsprep',0,irpnt,nr,ier)
  if (ier < 0) return
  call needw ('rsprep',1,iipnt,nr,ier)
  if (ier < 0) return
  maxlim = length/nr
  isym = 0
  if (nstore == 3) isym = 1
  call rsmatd (ndim,nr,nb,iwksp(maxnew),iwksp(jcnew),coef, &
    coef(ndim+1),coef(nr+ndim+1),coef(nr+1),maxrs, &
    iwksp(ijcrs),wksp(irs),maxlim,isym,ier)
  if (ier < 0) then
     call ershow (ier,'rsprep')
     return
  end if
  irpnt = irpnt + nr*maxrs
  iipnt = iipnt + maxrs
!
!  form reduced system rhs.
!
 45   irsrhs = irpnt
  ip1 = irpnt + nr
  call needw ('rsprep',0,irpnt,n+nr,ier)
  if (ier < 0) return
  if (nstore == 1) call rsbegp (n,nr,ndim,maxnz,jcoef, &
     coef,wksp(irsrhs),rhs,wksp(ip1))

  if (nstore >= 2) call rsrhsd (n,nr,ndim,iwksp(maxnew), &
    iwksp(jcnew),coef,wksp(irsrhs),rhs,wksp(ip1))

  irpnt = irpnt + nr
!
!  update constants.
!
  ndimrs = ndim
  mdimrs = mdim
  maxzrs = maxnz
  ndim = nr
  mdim = maxrs
  maxnz = maxrs
  nrr = nr
  return
end
subroutine rsrhsd (nn,nrr,ndim,maxnew,jcnew,coef,c,b,wksp)
!
!*******************************************************************************
!
!! RSRHSD computes  cr = br - t*inv(db)*bb.
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     diagonal storage
!
!  Parameters:
!
!        n          order of system
!                    systems)
!        nr         order of the red subsystem
!        ndim       row dimension of coef array
!        maxnew     number of columns in coef array
!        coef       real data structure
!        b          vector of length n containing bb behind br
!        c          vector of length nr containing cr
!        wksp       workspace vector of length nb
!
!  
!
  integer   jcnew(2,2), maxnew(2)
  dimension coef(ndim,2), b(1), c(1), wksp(1)
!
  n = nn
  nr = nrr
  nb = n - nr
!
!  compute cr.
!
  do 10 i = 1,nr
 10   c(i) = b(i)
  do 15 i = 1,nb
 15   wksp(i) = b(nr+i)/coef(nr+i,1)
  maxm1 = maxnew(1) - 1
  call vsubd (ndim,2,nr,nb,maxm1,coef(1,2),jcnew(1,2),c,wksp,nr)
  return
end
subroutine rsxbd (nn,nrr,ndim,maxnew,jcnew,coef,x,b)
!
!*******************************************************************************
!
!! RSXBD computes  xb = inv(db)*(bb - b*xr).
!
!                  a  = ( dr   t )
!                       ( b   db )
!
!     diagonal storage
!
!  Parameters:
!
!        n          order of system
!                    systems)
!        nr         order of the red subsystem
!        ndim       row dimension of coef array
!        maxnew     number of columns in coef array
!        coef       real data structure
!        x          vector of length n containing  xr, xb
!        b          vector of length n containing bb in the last
!                    nb locations
!
!  
!
  integer   jcnew(2,2), maxnew(2)
  dimension coef(ndim,2), x(1), b(1)
!
  n = nn
  nr = nrr
  nrp1 = nr + 1
  nb = n - nr
!
!  compute xb.
!
  do 10 i = nrp1,n
 10   x(i) = b(i)
  max2 = maxnew(2) - 1
  call vsubd (ndim,2,nb,nr,max2,coef(nrp1,2),jcnew(2,2),x(nrp1),x,-nr)
  do 15 i = nrp1,n
 15   x(i) = x(i)/coef(i,1)
  return
end
subroutine sbbs (ldd,ldt,n,kblszz,nsize,lbhb,iblock,d,t,jt,x,omega)
!
!*******************************************************************************
!
!! SBBS does an block SSOR backward pass.
!
!
!     symmetric diagonal data structure, natural ordering.
!     block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         kblsz    block size
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         lbhb     number of blocks per block row
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         x        input/output vector of length n
!         omega    over-relaxation factor
!
!  
!
  integer   jt(1), iblock(3,1)
  dimension d(ldd,2), t(ldt,1), x(1)
!
  kblsz = kblszz
  l = n/kblsz
  nt = iblock(3,1) - 1
  do 35 k = l,1,-1
     ist = (k - 1)*kblsz + 1
     ied = k*kblsz
     if (k == l) go to 15
     jjlim = min (lbhb,l-k+2)
     do 10 jj = 3,jjlim
        jblk = iblock(1,jj)
        jst = iblock(2,jj)
        mjj = iblock(3,jj)
        inc = jblk*kblsz
        istf = ist + inc
        if (istf > n) go to 10
        call vsubd (ldt,1,kblsz,kblsz,mjj,t(ist,jst),jt(jst),x(ist),x(istf),inc)
 10      continue
 15      if (nt >= 1) go to 25
     do 20 i = ist,ied
 20      x(i) = omega*d(i,1)*x(i)
     go to 35
 25      call bdsol (ldd,kblsz,nsize,nt,0,d(ist,1),x(ist),x(ist),0)
     do 30 i = ist,ied
 30      x(i) = omega*x(i)
 35   continue
  return
end
subroutine sbbsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  omega,iunif,wksp)
!
!*******************************************************************************
!
!! SBBSN does an block SSOR backward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   unif
!
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do backward solution.
!
 10   lm1 = l - 1
  do 50 k = lm1,1,-1
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     do 25 i = 1,na
 25      wksp(i) = 0.0
     do 30 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol <= k) go to 30
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        call vaddd (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb), &
                  wksp,x(istb),inc)
 30      continue
     if (ndt + ndb >= 1) go to 40
     do 35 i = ist,ied
 35      x(i) = x(i) - omega*d(i,1)*wksp(i-ist+1)
     go to 50
 40      call bdsol (ldd,na,nsize,ndt,ndb,d(ist,1),wksp,wksp,1)
     do 45 i = ist,ied
 45      x(i) = x(i) - omega*wksp(i-ist+1)
 50   continue
  return
end
subroutine sbbsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  omega,iunif)
!
!*******************************************************************************
!
!! SBBSNT does an block SSOR transpose backward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         omega    ove-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), x(1)
  logical   unif
!
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do backward solution.
!
 10   do 50 k = l,1,-1
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     if (ndt + ndb >= 1) go to 30
     do 25 i = ist,ied
 25      x(i) = omega*d(i,1)*x(i)
     go to 40
 30      call bdsolt (ldd,na,nsize,ndt,ndb,d(ist,1),x(ist),x(ist))
     do 35 i = ist,ied
 35      x(i) = omega*x(i)
 40      do 45 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol >= k) go to 45
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        call vsubdt (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb), &
                      x(istb),x(ist),inc)
 45      continue
 50   continue
  return
end
subroutine sbfs (ldd,ldt,n,kblszz,nsize,lbhb,iblock,d,t,jt,x,omega,wksp)
!
!*******************************************************************************
!
!! SBFS does an block SSOR forward pass.
!
!
!     symmetric diagonal data structure, natural ordering.
!     block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         kblsz    block size
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         lbhb     number of blocks per block row
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         x        input/output vector of length n
!         omega    over-relaxation factor
!         wksp     real workspace vector
!
!  
!
  integer   jt(1), iblock(3,1)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
!
  kblsz = kblszz
  l = n/kblsz
  lm1 = l - 1
  nt = iblock(3,1) - 1
  do 35 k = 1,lm1
     ist = (k - 1)*kblsz + 1
     ied = k*kblsz
     if (nt >= 1) go to 15
     do 10 i = ist,ied
 10      wksp(i-ist+1) = omega*d(i,1)*x(i)
     go to 25
 15      call bdsol (ldd,kblsz,nsize,nt,0,d(ist,1), x(ist),wksp,0)
     do 20 i = 1,kblsz
 20      wksp(i) = omega*wksp(i)
 25      jjlim = min (lbhb,l-k+2)
     do 30 jj = 3,jjlim
        jblk = iblock(1,jj)
        jst = iblock(2,jj)
        mjj = iblock(3,jj)
        inc = jblk*kblsz
        istf = ist + inc
        if (istf > n) go to 30
        call vsubdt (ldt,1,kblsz,kblsz,mjj,t(ist,jst),jt(jst), &
          x(istf),wksp,inc)
 30      continue
 35   continue
  return
end
subroutine sbfsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  omega,iunif)
!
!*******************************************************************************
!
!! SBFSN does an block SSOR forward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), x(1)
  logical   unif
!
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do forward solution.
!
 10   do 45 k = 1,l
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     do 25 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol >= k) go to 25
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        call vsubd (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb), &
          x(ist),x(istb),inc)
 25      continue
     if (ndt + ndb >= 1) go to 35
     do 30 i = ist,ied
 30      x(i) = omega*d(i,1)*x(i)
     go to 45
 35      call bdsol (ldd,na,nsize,ndt,ndb,d(ist,1),x(ist),x(ist),1)
     do 40 i = ist,ied
 40      x(i) = omega*x(i)
 45   continue
  return
end
subroutine sbfsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x, &
  omega,iunif,wksp)
!
!*******************************************************************************
!
!! SBFSNT does an block SSOR transpose forward solve.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         x        input/output vector of length n
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,2), t(ldt,1), wksp(1), x(1)
  logical   unif
!
  unif = iunif == 1
!
  l = ncolor
  if (.not. unif) go to 10
  na = nci(1)
  nb = na
  jlim = lbhb(1)
  l = n/na
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  kk = 1
!
!  do forward solution.
!
 10   lm1 = l - 1
  do 50 k = 1,lm1
     if (unif) go to 15
     kk = k
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     go to 20
 15      ist = (k - 1)*na + 1
 20      ied = ist + na - 1
     if (ndt + ndb >= 1) go to 30
     do 25 i = ist,ied
 25      wksp(i-ist+1) = omega*d(i,1)*x(i)
     go to 40
 30      call bdsolt (ldd,na,nsize,ndt,ndb,d(ist,1),x(ist),wksp)
     do 35 i = 1,na
 35      wksp(i) = omega*wksp(i)
 40      do 45 j = 3,jlim
        jcol = k + iblock(1,kk,j)
        if (jcol <= k) go to 45
        jstb = iblock(2,kk,j)
        mb = iblock(3,kk,j)
        if (unif) inc = (jcol - k)*na
        if (.not. unif) inc = ipt(jcol) - ipt(k)
        if (.not. unif) nb = nci(jcol)
        istb = ist + inc
        call vsubdt (ldt,ncolor,na,nb,mb,t(ist,jstb),jt(kk,jstb), &
          x(istb),wksp,inc)
 45      continue
 50   continue
  return
end
subroutine sbsl (ldd,ldt,n,kblsz,nsize,lbhb,iblock,d,t,jt,y,x,omega,wksp)
!
!*******************************************************************************
!
!! SBSL does an block SSOR solution.
!
!
!     symmetric diagonal data structure, natural ordering.
!     block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         kblsz    block size
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         lbhb     number of blocks per block row
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer vector giving the diagonal numbers
!                   for the off-diagonal blocks
!         y        input vector for the right-hand-side
!         x        output vector for the solution to q*x = y
!         omega    over-relaxation factor
!         wksp     real workspace vector
!
!  
!
  integer   jt(1), iblock(3,1)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  const = 2.0 - omega
  do 10 i = 1,n
 10   x(i) = const*y(i)
  call sbfs (ldd,ldt,n,kblsz,nsize,lbhb,iblock,d,t, jt,x,omega,wksp)
  call sbbs (ldd,ldt,n,kblsz,nsize,lbhb,iblock,d,t,jt,x,omega)
  return
end
subroutine sbsln (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb, &
  iblock,d,t,jt,y,x,omega,iunif,wksp)
!
!*******************************************************************************
!
!! SBSLN does an block SSOR solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  const = 2.0 - omega
  do 10 i = 1,n
 10   x(i) = const*y(i)
  call sbfsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega,iunif)
  call sbbsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega, &
     iunif,wksp)
  return
end
subroutine sbsln1 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  omega,iunif)
!
!*******************************************************************************
!
!! SBSLN1 does an block SSOR forward solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1),iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), x(1), y(1)
!
  const = 2.0 - omega
  do 10 i = 1,n
 10   x(i) = const*y(i)
  call sbfsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega,iunif)
  return
end
subroutine sbsln2 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  omega,iunif,wksp)
!
!*******************************************************************************
!
!! SBSLN2 does an block SSOR back solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1),iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call sbbsn (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega, &
    iunif,wksp)
  return
end
subroutine sbsln3 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  omega,iunif)
!
!*******************************************************************************
!
!! SBSLN3 does an block SSOR transpose forward solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), x(1), y(1)
!
  const = 2.0 - omega
  do 10 i = 1,n
 10   x(i) = const*y(i)
  call sbbsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega,iunif)
  return
end
subroutine sbsln4 (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  omega,iunif,wksp)
!
!*******************************************************************************
!
!! SBSLN4 does an block SSOR transpose back solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1),iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call sbfsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega, &
    iunif,wksp)
  return
end
subroutine sbslnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,y,x, &
  omega,iunif,wksp)
!
!*******************************************************************************
!
!! SBSLNT does an block SSOR transpose solution.
!
!
!     nonsymmetric diagonal data structure, natural or multi-color
!     orderings, block ssor preconditioning.
!
!  Parameters:
!
!         ldd      row dimension of d array
!         ldt      row dimension of t array
!         n        size of system
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!                   ncolor = 1 if iunif = 1.
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!                   if iunif = 1, nci(1) is the constant block size.
!         ipt      integer pointer vector of length ncolor+1 if
!                   iunif = 0.  formed in the factorization routine.
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!                   if iunif = 1, lbhb is of length 1.
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         d        array for diagonal block
!         t        array for off-diagonal blocks
!         jt       integer array of size ncolor by whatever
!                   giving the off-diagonal block diagonal numbers
!                   for each distinct block size.  jd is 1 by whatever
!                   if iunif = 1.
!         y        input vector of length n containing right-hand-side
!         x        output vector containing the solution to q*x = y
!         omega    over-relaxation factor
!         iunif    uniform block size switch
!                   = 0   diagonal blocks are not of uniform size
!                   = 1   diagonal blocks are of uniform size
!         wksp     real workspace vector
!
!  
!
  integer   ipt(1), jt(ncolor,1), nci(1), lbhb(1), iblock(3,ncolor,2)
  dimension d(ldd,1), t(ldt,1), wksp(1), x(1), y(1)
!
  const = 2.0 - omega
  do 10 i = 1,n
 10   x(i) = const*y(i)
  call sbfsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega, &
    iunif,wksp)
  call sbbsnt (ldd,ldt,n,nsize,ncolor,nci,ipt,lbhb,iblock,d,t,jt,x,omega,iunif)
  return
end
subroutine scal1 (nn,ndim,maxnzz,jcoef,coef,rhs,u,ubar,diag,work,iflag,ier)
!
!*******************************************************************************
!
!! SCAL1 scales the original matrix to a unit diagonal matrix.
!
!
!     (Purdue data structure)
!     rhs and u vectors are scaled accordingly.  upon output, diag
!     contains the reciprocal square roots of the diagonal elements.
!     it is assumed that the diagonal of the matrix is in column one
!     of coef.
!
!  Parameters:
!
!         n       dimension of matrix
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         ubar    exact solution (optional)
!         diag    work array of length n (nonvolatile)
!         work    work array of length n (volatile)
!         iflag   flag for ubar
!                  = 0  do not scale ubar
!                  = 1  scale ubar
!         ier     error flag -- on return, values mean
!                      0 -- no errors detected
!                     -4 -- nonpositive diagonal element
!
!  
!
  integer   jcoef(ndim,1)
  dimension coef(ndim,1), rhs(1), u(1), diag(1), work(1), ubar(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
  maxnz = maxnzz
!
!  check for positive diagonal entries for each row.
!
  cmin = vmin (n,coef)
  if (cmin > 0.0) go to 10
!
!  fatal error -- nonpositive diagonal element.
!
  ier = -4
  return
!
!  scale matrix.  store reciprocal square roots
!  of diagonal entries in diag.
!
 10   do 15 i = 1,n
 15   diag(i) = sqrt (coef(i,1))
!
!  scale rhs, u, and ubar.
!
  do 20 i = 1,n
 20   u(i) = diag(i)*u(i)
  if (iflag == 0) go to 30
  do 25 i = 1,n
 25   ubar(i) = diag(i)*ubar(i)
 30   do 35 i = 1,n
 35   diag(i) = 1.0/diag(i)
  do 40 i = 1,n
 40   rhs(i) = diag(i)*rhs(i)
  if (keygs == 2) go to 55
!
!  using gathers.
!
  do 50 j = 1,maxnz
     call vgathr (n,diag,jcoef(1,j),work)
     do 45 i = 1,n
 45      coef(i,j) = diag(i)*coef(i,j)*work(i)
 50   continue
  return
!
!  not using gathers.
!
 55   do 65 j = 1,maxnz
     do 60 i = 1,n
 60      coef(i,j) = diag(i)*coef(i,j)*diag(jcoef(i,j))
 65   continue
  return
end
subroutine scal2 (nn,ndim,maxnz,jcoef,coef,rhs,u,ubar,diag,iflag,ier)
!
!*******************************************************************************
!
!! SCAL2 scales the original matrix to a unit diagonal matrix.
!
!
!     (diagonal data structure)
!     rhs and u vectors are scaled accordingly.  upon output, diag
!     contains the reciprocal square roots of the diagonal elements.
!     it is assumed that the diagonal of the matrix is in column one
!     of coef.
!
!  Parameters:
!
!         n       dimension of matrix
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         ubar    exact solution (optional)
!         diag    work array of length n (nonvolatile)
!         iflag   flag for ubar
!                  = 0  do not scale ubar
!                  = 1  scale ubar
!         ier     error flag -- on return, values mean
!                      0 -- no errors detected
!                     -4 -- nonpositive diagonal element
!
!  
!
  integer   jcoef(2)
  dimension coef(ndim,1), rhs(1), u(1), diag(1), ubar(1)
!
!
  n = nn
!
!  check for positive diagonal entries for each row.
!
  cmin = vmin (n,coef)
  if (cmin > 0.0) go to 10
!
!  fatal error -- nonpositive diagonal element.
!
  ier = -4
  return
!
!  scale matrix.  store reciprocal square roots
!  of diagonal entries in diag.
!
 10   do 15 i = 1,n
 15   diag(i) = sqrt (coef(i,1))
!
!  scale rhs, u, and ubar.
!
  do 20 i = 1,n
 20   u(i) = diag(i)*u(i)
  if (iflag == 0) go to 30
  do 25 i = 1,n
 25   ubar(i) = diag(i)*ubar(i)
 30   do 35 i = 1,n
 35   diag(i) = 1.0/diag(i)
  do 40 i = 1,n
 40   rhs(i) = diag(i)*rhs(i)
!
!  scale matrix.
!
  do 60 j = 1,maxnz
     ind = jcoef(j)
     len = n - iabs(ind)
     if (ind < 0) go to 50
     do 45 i = 1,len
 45      coef(i,j) = diag(i)*coef(i,j)*diag(i+ind)
     go to 60
 50      do 55 i = 1,len
 55      coef(i-ind,j) = diag(i-ind)*coef(i-ind,j)*diag(i)
 60   continue
  return
end
subroutine scal3 (nn,nz,ia,ja,a,rhs,u,ubar,diag,work,iflag,ier)
!
!*******************************************************************************
!
!! SCAL3 scales the original matrix to a unit diagonal matrix.
!
!
!     (sparse data structure)
!     rhs and u vectors are scaled accordingly.  upon output, diag
!     contains the reciprocal square roots of the diagonal elements.
!     it is assumed that the diagonal of the matrix is in the
!     n first locations of a.
!
!  Parameters:
!
!         n       dimension of matrix
!         nz      length of ia, ja, and a vectors
!         a       vector containing matrix coefficients
!         ia      vector of i values
!         ja      vector of j values
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         ubar    exact solution (optional)
!         diag    vector of length n containing the reciprocal
!                  square roots of the diagonal elements upon
!                  output
!         work    workspace vector of length n
!         iflag   flag for ubar
!                  = 0  do not scale ubar
!                  = 1  scale ubar
!         ier     error flag -- on return, values mean
!                      0 -- no errors detected
!                     -4 -- nonpositive diagonal element
!
!  
!
  integer   ia(1), ja(1)
  dimension a(1), rhs(1), u(1), diag(1), work(1), ubar(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
!
!  check for positive diagonal entries for each row.
!
  cmin = vmin (n,a)
  if (cmin > 0.0) go to 10
!
!  fatal error -- nonpositive diagonal element.
!
  ier = -4
  return
!
!  scale matrix.  store reciprocal square roots
!  of diagonal entries in diag.
!
 10   do 15 i = 1,n
 15   diag(i) = sqrt (a(i))
!
!  scale rhs, u, and ubar.
!
  do 20 i = 1,n
 20   u(i) = diag(i)*u(i)
  if (iflag == 0) go to 30
  do 25 i = 1,n
 25   ubar(i) = diag(i)*ubar(i)
 30   do 35 i = 1,n
 35   diag(i) = 1.0/diag(i)
  do 40 i = 1,n
 40   rhs(i) = diag(i)*rhs(i)
  if (keygs == 2) go to 60
!
!  using gathers.
!
  ist = 1
 45   ied = min (ist-1+n,nz)
  if (ied < ist) return
     len = ied - ist + 1
     call vgathr (len,diag,ia(ist),work)
     do 50 i = ist,ied
 50      a(i) = a(i)*work(i-ist+1)
     call vgathr (len,diag,ja(ist),work)
     do 55 i = ist,ied
 55      a(i) = a(i)*work(i-ist+1)
  ist = ied + 1
  go to 45
!
!  not using gathers.
!
 60   do 65 i = 1,nz
 65   a(i) = a(i)*diag(ia(i))*diag(ja(i))
  return
end
subroutine scale (coef,jcoef,wksp,icall,n,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! SCALE scales the matrix, U, UBAR, and RHS.
!
!  Parameters:
!
!       icall     key to indicate whether scaling (icall=1) or
!                  unscaling (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -4   nonpositive diagonal element
!
!  
!
  integer jcoef(2)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!  data common blocks
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  if (iscale /= 1) return
  go to (5,10,10,15,15), nstore
 5    call scalep (coef,jcoef,wksp,icall,n,u,ubar,rhs,ier)
  return
 10   call scaled (coef,jcoef,wksp,icall,n,u,ubar,rhs,ier)
  return
 15   call scales (coef,jcoef,wksp,icall,n,u,ubar,rhs,ier)
  return
end
subroutine scaled (coef,jcoef,wksp,icall,nn,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! SCALED scales the matrix, U, UBAR, and RHS.
!
!
!     (symmetric or nonsymmetric diagonal format)
!
!  Parameters:
!
!       icall     key to indicate whether scaling (icall=1) or
!                  unscaling (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -4   nonpositive diagonal element
!
!  
!
  integer jcoef(2)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  iflag = 0
  if (ntest == 6) iflag = 1
  if (icall == 2) go to 20
!
!  scale system.
!
!  check for sufficient room.
!
  call needw ('scaled',0,irpnt,n,ier)
  if (ier < 0) return
  iptscl = irpnt
  irpnt = irpnt + n
  call scal2 (n,ndim,maxnz,jcoef,coef,rhs,u,ubar,wksp(iptscl),iflag,ier)
  if (ier < 0) call ershow (ier,'scaled')
  return
!
!  unscale system.
!
 20   call uscal2 (n,ndim,maxnz,jcoef,coef,rhs,u,ubar,wksp(iptscl),iflag)
  return
end
subroutine scalep (coef,jcoef,wksp,icall,nn,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! SCALEP scales the matrix, U, UBAR, and RHS.
!
!
!     (Purdue format)
!
!  Parameters:
!
!       icall     key to indicate whether scaling (icall=1) or
!                  unscaling (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -4   nonpositive diagonal element
!
!  
!
  integer jcoef(2)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  iflag = 0
  if (ntest == 6) iflag = 1
  if (icall == 2) go to 20
!
!  scale system.
!
!  check for sufficient room.
!
  call needw ('scalep',0,irpnt,2*n,ier)
  if (ier < 0) return
  iptscl = irpnt
  irpnt = irpnt + n
  call scal1 (n,ndim,maxnz,jcoef,coef,rhs,u,ubar,wksp(iptscl),wksp(irpnt), &
    iflag,ier)
  if (ier < 0) call ershow (ier,'scalep')
  return
!
!  unscale system.
!
 20   call uscal1 (n,ndim,maxnz,jcoef,coef,rhs,u,ubar,wksp(iptscl), &
    wksp(irpnt),iflag)
  return
end
subroutine scales (coef,jcoef,wksp,icall,nn,u,ubar,rhs,ier)
!
!*******************************************************************************
!
!! SCALES scales the matrix, U, UBAR, and RHS.  (sparse format)
!
!  Parameters:
!
!       icall     key to indicate whether scaling (icall=1) or
!                  unscaling (icall=2) is to be done
!       n         order of system
!       u         current solution estimate
!       ubar      input vector containing the true solution
!                  (optional)
!       rhs       input vector.  contains the right hand side
!                 of the matrix problem.
!       ier       error flag
!                  =   0   no errors detected
!                  =  -4   nonpositive diagonal element
!
!  
!
  integer jcoef(2)
  dimension rhs(1), u(1), ubar(1), coef(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
!
!
!
!  data common blocks
!
  common / dscons / ndim, mdim, maxnz
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  iflag = 0
  if (ntest == 6) iflag = 1
  if (icall == 2) go to 10
!
!  scale system.
!
!  check for sufficient room.
!
  call needw ('scales',0,irpnt,2*n,ier)
  if (ier < 0) return
  iptscl = irpnt
  irpnt = irpnt + n
  call scal3 (n,maxnz,jcoef,jcoef(ndim+1),coef,rhs,u,ubar, &
    wksp(iptscl),wksp(irpnt),iflag,ier)
  if (ier < 0) call ershow (ier,'scales')
  return
!
!  unscale system.
!
 10   call uscal3 (n,maxnz,jcoef,jcoef(ndim+1),coef,rhs,u,ubar,wksp(iptscl), &
    wksp(irpnt),iflag)
  return
end
function seldot (n,iu,u1,u2,u3,u4,iv,v1,v2,v3,v4)
!
!*******************************************************************************
!
!! SELDOT computes a dot product from a selected pair of vectors.
!
  dimension u1(1), u2(1), u3(1), v1(1), v2(1), v3(1)
  dimension u4(1), v4(1)
!
  ind = 1 + (iv-1) + 4*(iu-1)
  if (ind == 1)  seldot = vdot (n,u1,v1)
  if (ind == 2)  seldot = vdot (n,u1,v2)
  if (ind == 3)  seldot = vdot (n,u1,v3)
  if (ind == 4)  seldot = vdot (n,u1,v4)
  if (ind == 5)  seldot = vdot (n,u2,v1)
  if (ind == 6)  seldot = vdot (n,u2,v2)
  if (ind == 7)  seldot = vdot (n,u2,v3)
  if (ind == 8)  seldot = vdot (n,u2,v4)
  if (ind == 9)  seldot = vdot (n,u3,v1)
  if (ind == 10) seldot = vdot (n,u3,v2)
  if (ind == 11) seldot = vdot (n,u3,v3)
  if (ind == 12) seldot = vdot (n,u3,v4)
  if (ind == 13) seldot = vdot (n,u4,v1)
  if (ind == 14) seldot = vdot (n,u4,v2)
  if (ind == 15) seldot = vdot (n,u4,v3)
  if (ind == 16) seldot = vdot (n,u4,v4)
  return
end
subroutine si (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SI is the user interface to the Chebyshev acceleration algorithm.
!
  external suba, subat, subql, subqlt, subqr, subqrt, subadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  ier = 0
  call needw ('si',0,irpnt,4*n,ier)
  if (ier < 0) return
  nw = lenr - irpnt + 1
  call siw (suba,subql,coef,jcoef,wksp,iwksp,n,u,ubar,rhs,wksp(irpnt),nw, &
    iparm,rparm,ier)
  irmax = irpnt + nw - 1
  return
end
subroutine siw (suba,subq,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs,wksp,nw,iparm, &
  rparm,ier)
!
!*******************************************************************************
!
!! SIW drives the Chebyshev acceleration algorithm.
!
!  Parameters:
!
!          suba   matrix-vector multiplication routine
!          subq   preconditioning routine
!          n      input integer.  order of the system (= nn)
!          u      input/output vector.  on input, u contains the
!                 initial guess to the solution.  on output, it
!                 contains the latest estimate to the solution.
!          ubar   input vector containing the true solution
!                  (optional)
!          rhs    input vector.  contains the right hand side
!                 of the matrix problem.
!          wksp   vector used for working space.
!          nw     length of wksp array.  if this length is less than
!                  the amount needed, nw will give the needed amount
!                  upon output.
!          iparm  integer vector of length 30.  allows user to
!                 specify some integer parameters which affect
!                 the method.
!          rparm  real vector of length 30.  allows user to
!                 specify some real parameters which affect
!                 the method.
!          ier    output integer.  error flag.
!
!  
!
  external  suba, subq
  integer   iparm(30), jcoef(2), jwfac(1)
  dimension rhs(1), u(1), ubar(1), wksp(1), rparm(30), coef(1),wfac(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
!  initialize common blocks
!
  ier = 0
  n = nn
  t1 = timer (dummy)
  iacel = 2
  timit = 0.0
  digit1 = 0.0
  digit2 = 0.0
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 35
  if (level >= 2) write (nout,10)
 10   format (1x,'si')
!
!  compute workspace base addresses and check for sufficient
!  workspace.
!
  iw1 = 1
  iw2 = iw1 + n
  iw3 = iw2 + n
  iw4 = iw3 + n
  nwksp = 4*n
  if (nw >= nwksp) go to 15
  ier = -2
  call ershow (ier,'siw')
  go to 30
 15   continue
  call nmcalc (coef,jcoef,wfac,jwfac,1,subq,n,rhs,ubar,wksp,ier)
  if (ier < 0) go to 30
!
!  compute an initial rayleigh quotient and adjust emax, emin.
!
  call vfill (n,wksp,1.0)
  call subq (coef,jcoef,wfac,jwfac,n,wksp,wksp(iw2))
  call suba (coef,jcoef,wfac,jwfac,n,wksp(iw2),wksp(iw3))
  rq = vdot (n,wksp(iw2),wksp(iw3)) /vdot (n,wksp(iw2),wksp)
  rqmax = rq
  rqmin = rq
  if (maxadd) emax = amax1 (emax,rqmax)
  if (minadd) emin = amin1 (emin,rqmin)
  if (minadd) emin = amax1 (emin,0.0)
!
!  zero out workspace
!
  call vfill (nwksp,wksp,0.0)
!
!  iteration sequence
!
  call itsi (suba,subq,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    wksp(iw1),wksp(iw2),wksp(iw3),wksp(iw4),ier)
!
  if (ier < 0  .or.  ier == 1) go to 25
!
!  method has converged
!
  if (level >= 1) write (nout,20) in
 20   format (/1x,'si  has converged in ',i5,' iterations ')
!
!  optional error analysis
!
 25   if (idgts < 0) go to 30
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wksp,digit1,digit2,idgts)
!
!  set return parameters in iparm and rparm
!
 30   t2 = timer (dummy)
  nw = 4*n
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
  rparm(9) = omega
  rparm(10) = alphab
  rparm(11) = betab
  rparm(12) = specr
!
 35   continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
!
  return
end
subroutine sor (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SOR is the user interface to the SOR algorithm.
!
  external suba, subat, subql, subqlt, subqr, subqrt, subadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  ier = 0
  call needw ('sor',0,irpnt,2*n,ier)
  if (ier < 0) return
  nw = lenr - irpnt + 1
  call sorw (suba,subql,coef,jcoef,wksp,iwksp, &
    n,u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = irpnt + nw - 1
  return
end
subroutine sor1 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SOR1 drives the point SOR method.
!
  external accel, suba8, suba9, subq78, noadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call move1 (ndim,mdim,n,maxnz,jcoef,coef,maxt,maxb,ier)
  if (ier < 0) then
     call ershow (ier,'sor1')
     return
  end if
  call split (accel,suba8,suba9,subq78,subq78,subq78,subq78,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine sor2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SOR2 drives the point SOR method.
!
  external accel, suba1, subq6, noadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call rowise (maxnz,jcoef,irwise)
  call needw ('sor2',1,iipnt,maxnz,ier)
  if (ier < 0) return
  iwkpt1 = iipnt
  iipnt = iipnt + maxnz
  call split (accel,suba1,suba1,subq6,subq6,subq6,subq6,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  iipnt = iipnt - maxnz
  return
end
subroutine sor3 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SOR3 drives the point SOR method.
!
  external accel, suba4, suba5, subq40, noadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call rowise (maxnz,jcoef,irwise)
  call needw ('sor3',1,iipnt,maxnz,ier)
  if (ier < 0) return
  call needw ('sor3',0,irpnt,n,ier)
  if (ier < 0) return
  call move2 (ndim,n,maxnz,jcoef,coef,wksp(irpnt),iwksp(iipnt),maxt,maxb)
  iwkpt1 = iipnt
  iipnt = iipnt + maxnz
  call split (accel,suba4,suba5,subq40,subq40,subq40,subq40,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  iipnt = iipnt - maxnz
  return
end
subroutine sor6 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SOR6 drives the multi-color SOR method.
!
  external accel, suba8, subq96, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba8,suba8,subq96,subq96,subq96,subq96,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine sor7 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SOR7 drives the multi-color SOR method.
!
  external accel, suba2, subq26, copy, noadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
!
  t1 = timer (dummy)
  if (ifact == 1) call mfact (coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  call split (accel,suba2,suba2,subq26,subq26,subq26,subq26,copy,copy,noadp, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  return
end
subroutine sorcp (ndimm,n,jc,d,c,ncol,nc,nt,nb,omega,u,rhs,unew)
!
!*******************************************************************************
!
!! SORCP does an SOR solve.
!     (Purdue storage, multicolor)
!
!        unew = inv((1/w)*d + l)*(((1-w)/w)*d*un + (rhs - u*un))
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          omega  over-relaxation factor
!          u      current solution
!          rhs    right-hand-side
!          unew   updated solution
!
!  
!
  integer   jc(ndimm,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimm,1), u(1), rhs(1), unew(1)
!
  ndim = ndimm
  ncolor = ncol
!
!  rhs = ((1-w)/w)*d*un + (rhs - u*un)
!
  ist =  1
  do 10 icol = 1,ncolor
     npt = nc(icol)
     j2 = nt(icol)
     call vsubp (ndim,ndim,npt,j2,c(ist,1),jc(ist,1),rhs(ist),u,unew)
     ist = ist + npt
 10   continue
  con = (1.0 - omega)/omega
  do 15 i = 1,n
 15   unew(i) = con*d(i)*u(i) + rhs(i)
!
!  unew = inv((1/w)*d + l)*rhs
!
  ist = 1
  do 25 icol = 1,ncolor
     npt = nc(icol)
     ied = ist + npt - 1
     j1 = nt(icol) + 1
     mj = nb(icol)
     call vsubp (ndim,ndim,npt,mj,c(ist,j1),jc(ist,j1),unew(ist),unew,rhs)
     do 20 i = ist,ied
 20      unew(i) = omega*unew(i)/d(i)
     ist = ist + npt
 25   continue
  return
end
subroutine sordb (ldf,ndim,nsize,kblszz,iblock,lbhb,dfac,coef,jcoef,nn, &
  omega,u,rhs,unew)
!
!*******************************************************************************
!
!! SORDB does an SOR pass
!     (symmetric block diagonal format, constant block size)
!
!        unew = inv((1/w)*d + l)*(((1-w)/w)*d*un + (rhs - u*un))
!
!  Parameters:
!
!         ldf      row dimension of dfac
!         ndim     row dimension of coef array
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         lbhb     column size of iblock
!         dfac     array for diagonal block factorization
!         coef     array for matrix coefficients
!         jcoef    vector for diagonal numbers
!         n        size of system
!         omega    relaxation parameter
!         u        current solution estimate
!         rhs      right-hand-side
!         unew     updated solution estimate
!
!  
!
  integer   jcoef(2), iblock(3,1)
  dimension dfac(ldf,1), coef(ndim,2), u(1), rhs(1), unew(1)
!
  n = nn
  kblsz = kblszz
!
!  rhs = ((1-w)/w)*d*un + (rhs - u*un)
!
  nwdiag = iblock (3,1)
  nt = nwdiag - 1
  maxt = 0
  if (lbhb < 3) go to 15
  do 10 j = 3,lbhb
     maxt = maxt + iblock(3,j)
 10   continue
 15   jbgn = nwdiag + 1
  call vsubd (ndim,1,n,n,maxt,coef(1,jbgn),jcoef(jbgn),rhs,u,0)
  call bmul (ndim,n,nt,coef,coef(1,2),u,unew)
  con = (1.0 - omega)/omega
  do 20 i = 1,n
 20   unew(i) = con*unew(i) + rhs(i)
!
!  unew = inv((1/w)*d + l)*rhs
!
  l = n/kblsz
  do 50 k = 1,l
     ist = (k - 1)*kblsz + 1
     ied = k*kblsz
     if (nt >= 1) go to 30
     do 25 i = ist,ied
 25      unew(i) = omega*dfac(i,1)*unew(i)
     go to 40
 30      call bdsol (ldf,kblsz,nsize,nt,0,dfac(ist,1),unew(ist),unew(ist),0)
     do 35 i = ist,ied
 35      unew(i) = omega*unew(i)
 40      if (k == l) go to 50
     jjlim = min (lbhb,l-k+2)
     do 45 jj = 3,jjlim
        jblk = iblock(1,jj)
        jst = iblock(2,jj) + nwdiag
        mjj = iblock(3,jj)
        inc = jblk*kblsz
        istf = ist + inc
        if (istf > n) go to 45
        call vsubdt (ndim,1,kblsz,kblsz,mjj,coef(ist,jst), &
          jcoef(jst),unew(istf),unew(ist),inc)
 45      continue
 50   continue
  return
end
subroutine sordmb (ldf,ndim,nsize,iblock,lbhb,ncol,nc,ipt,dfac,coef,jcnew,nn, &
  omega,u,rhs,unew)
!
!*******************************************************************************
!
!! SORDMB does an SOR pass.
!     (nonsymmetric block diagonal format, nonconstant block size)
!
!        unew = inv((1/w)*d + l)*(((1-w)/w)*d*un + (rhs - u*un))
!
!  Parameters:
!
!         ldf      row dimension of dfac array
!         ndim     row dimension of coef array
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!         ncolor   number of distinct block sizes
!         nc       integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!         ipt      integer pointer vector of length ncolor+1
!                   giving the starting locations of new block
!                   rows
!         dfac     array for diagonal block factorization
!         coef     array of matrix coefficients
!         jcnew    integer array of row dimension ncolor giving the
!                   diagonal numbers for each block
!         n        size of system
!         omega    relaxation parameter
!         u        current solution estimate
!         rhs      right-hand-side
!         unew     updated solution estimate
!
!  
!
  integer   jcnew(ncol,1), iblock(3,ncol,2), lbhb(1), nc(1),ipt(1)
  dimension dfac(ldf,1), coef(ndim,2), u(1), rhs(1), unew(1)
!
  n = nn
  ncolor = ncol
!
!  rhs = ((1-w)/w)*d*un + (rhs - u*un)
!
  ndt = iblock (3,1,1) - 1
  ndb = iblock (3,1,2)
  nwdiag = ndt + ndb + 1
  do 15 k = 1,ncolor
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nc(k)
     do 10 j = 3,jlim
        jcol = k + iblock(1,k,j)
        if (jcol <= k .or. jcol > ncolor) go to 10
        jstb = iblock(2,k,j) + nwdiag
        mb = iblock(3,k,j)
        inc = ipt(jcol) - ipt(k)
        nb = nc(jcol)
        istb = ist + inc
        call vsubd (ndim,ncolor,na,nb,mb,coef(ist,jstb),jcnew(k,jstb), &
          rhs(ist),u(istb),inc)
 10      continue
 15   continue
  ind = ndt + 2
  call bmuln (ndim,n,ndt,ndb,coef,coef(1,2),coef(1,ind),u,unew)
  con = (1.0 - omega)/omega
  do 20 i = 1,n
 20   unew(i) = con*unew(i) + rhs(i)
!
!  unew = inv((1/w)*d + l)*rhs
!
  do 45 k = 1,ncolor
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nc(k)
     ndt = iblock(3,k,1) - 1
     ndb = iblock(3,k,2)
     ied = ist + na - 1
     do 25 j = 3,jlim
        jcol = k + iblock(1,k,j)
        if (jcol >= k .or. jcol <= 0) go to 25
        jstb = iblock(2,k,j) + nwdiag
        mb = iblock(3,k,j)
        inc = ipt(jcol) - ipt(k)
        nb = nc(jcol)
        istb = ist + inc
        call vsubd (ndim,ncolor,na,nb,mb,coef(ist,jstb), &
          jcnew(k,jstb),unew(ist),unew(istb),inc)
 25      continue
     if (ndt + ndb >= 1) go to 35
     do 30 i = ist,ied
 30      unew(i) = omega*dfac(i,1)*unew(i)
     go to 45
 35      call bdsol (ldf,na,nsize,ndt,ndb,dfac(ist,1),unew(ist),unew(ist),1)
     do 40 i = ist,ied
 40      unew(i) = omega*unew(i)
 45   continue
  return
end
subroutine sordn (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,omegaa,irwise,u,rhs,unew, &
  iwksp)
!
!*******************************************************************************
!
!! SORDN does an SOR solve (natural ordering, nonsymmetric diagonal storage).
!
!        unew = inv(d + w*l)*((1-w)*d*un + w*(rhs - u*un))
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        b      array of active size n by maxb giving the sub-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        u      current solution vector
!        rhs    right hand side
!        unew   updated solution vector
!        iwksp  integer workspace of length maxt
!
!  
!
  dimension d(1), t(ndim,1), b(ndim,1), u(1), unew(1), rhs(1)
  integer   jt(1), jb(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  omega = omegaa
!
!  rhs = (1-w)*d*un + w*(rhs - u*un)
!
  call vsubd (ndim,1,n,n,maxt,t,jt,rhs,u,0)
  con = 1.0 - omega
  do 10 i = 1,n
 10   rhs(i) = con*d(i)*u(i) + omega*rhs(i)
!
!  rhs = inv(i+w*l*inv(d))*rhs
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 50
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxb
 15   iwksp(i) = 1 - jb(i)
!
!  determine nc, imin.
!
 20   nc = n
  do 25 i = 1,maxb
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 25
     nc = nterm
     imin = i
 25   continue
  if (nc >= n) go to 70
  ndel = -jb(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 40
!
!  special case for first minor subdiagonal.
!
  nc1 = n
  do 30 i = 1,maxb
     if (i == imin) go to 30
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imin) = nc1 + 1
  do 35 j = ibeg,nc1
 35   rhs(j) = rhs(j) - omega*b(j,imin)*rhs(j-1)/d(j-1)
  go to 20
!
!  far diagonals  (do vector computations).
!
 40   iwksp(imin) = iwksp(imin) + ndel
  iend = min (ibeg+ndel-1,n)
  do 45 i = ibeg,iend
 45   rhs(i) = rhs(i) - omega*b(i,imin)*rhs(i-ndel)/d(i-ndel)
  go to 20
!
!  rowwise algorithm.
!
 50   do 65 i = 1,n
     do 55 j = 1,maxb
 55      iwksp(j) = max (1,i+jb(j))
     sum = 0.0
     do 60 j = 1,maxb
 60      sum = sum + b(i,j)*rhs(iwksp(j))/d(iwksp(j))
     rhs(i) = rhs(i) - omega*sum
 65   continue
!
!  unew = inv(d)*rhs
!
 70   do 75 i = 1,n
 75   unew(i) = rhs(i)/d(i)
  return
end
subroutine sordnb (ldf,ndim,nsize,kblszz,iblock,lbhbb,dfac,coef,jcoef,nn, &
  omega,u,rhs,unew)
!
!*******************************************************************************
!
!! SORDNB does an SOR pass.
!     (nonsymmetric block diagonal format, constant block size)
!
!        unew = inv((1/w)*d + l)*(((1-w)/w)*d*un + (rhs - u*un))
!
!  Parameters:
!
!         ldf      row dimension of dfac
!         ndim     row dimension of coef array
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         lbhb     column size of iblock
!         dfac     array for diagonal block factorization
!         coef     array for matrix coefficients
!         jcoef    vector for diagonal numbers
!         n        size of system
!         omega    relaxation parameter
!         u        current solution estimate
!         rhs      right-hand-side
!         unew     updated solution estimate
!
!  
!
  integer   jcoef(2), iblock(3,2)
  dimension dfac(ldf,1), coef(ndim,2), u(1), rhs(1), unew(1)
!
  n = nn
  kblsz = kblszz
  lbhb = lbhbb
!
!  rhs = ((1-w)/w)*d*un + (rhs - u*un)
!
  nt = iblock (3,1) - 1
  nb = iblock (3,2)
  nwdiag = nt + nb + 1
  maxt = 0
  if (lbhb < 3) go to 15
  do 10 j = 3,lbhb
     ind = iblock(1,j)
     if (ind > 0) maxt = maxt + iblock(3,j)
 10   continue
 15   jbgn = nwdiag + 1
  call vsubd (ndim,1,n,n,maxt,coef(1,jbgn),jcoef(jbgn),rhs,u,0)
  ind = nt + 2
  call bmuln (ndim,n,nt,nb,coef,coef(1,2),coef(1,ind),u,unew)
  con = (1.0 - omega)/omega
  do 20 i = 1,n
 20   unew(i) = con*unew(i) + rhs(i)
!
!  unew = inv((1/w)*d + l)*rhs
!
  l = n/kblsz
  do 45 k = 1,l
     ist = (k - 1)*kblsz + 1
     ied = k*kblsz
     do 25 j = 3,lbhb
        jcol = k + iblock(1,j)
        if (jcol >= k .or. jcol <= 0) go to 25
        jstb = iblock(2,j) + nwdiag
        mb = iblock(3,j)
        inc = (jcol - k)*kblsz
        istb = ist + inc
        call vsubd (ndim,1,kblsz,kblsz,mb,coef(ist,jstb), &
          jcoef(jstb),unew(ist),unew(istb),inc)
 25      continue
     if (nt + nb >= 1) go to 35
     do 30 i = ist,ied
 30      unew(i) = omega*dfac(i,1)*unew(i)
     go to 45
 35      call bdsol (ldf,kblsz,nsize,nt,nb,dfac(ist,1),unew(ist),unew(ist),1)
     do 40 i = ist,ied
 40      unew(i) = omega*unew(i)
 45   continue
  return
end
subroutine sords (ndim,nn,maxtt,jt,d,t,omegaa,irwise,u,rhs,unew,iwksp)
!
!*******************************************************************************
!
!! SORDS does an SOR solve (natural ordering,
!     symmetric diagonal storage).
!
!        unew = inv(d + w*l)*((1-w)*d*un + w*(rhs - u*un))
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        u      current solution vector
!        rhs    right hand side
!        unew   updated solution vector
!        iwksp  integer workspace of length maxt
!
!  
!
  dimension d(1), t(ndim,1), u(1), unew(1), rhs(1)
  integer   jt(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  omega = omegaa
!
!  rhs = (1-w)*d*un + w*(rhs - u*un)
!
  call vsubd (ndim,1,n,n,maxt,t,jt,rhs,u,0)
  con = 1.0 - omega
  do 10 i = 1,n
 10   rhs(i) = con*d(i)*u(i) + omega*rhs(i)
!
!  rhs = inv(i+w*l*inv(d))*rhs
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 50
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxt
 15   iwksp(i) = jt(i) + 1
!
!  determine nc, imin.
!
 20   nc = n
  do 25 i = 1,maxt
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 25
     nc = nterm
     imin = i
 25   continue
  if (nc >= n) go to 70
  ndel = jt(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 40
!
!  special case for first minor subdiagonal.
!
  nc1 = n
  do 30 i = 1,maxt
     if (i == imin) go to 30
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imin) = nc1 + 1
  do 35 j = ibeg,nc1
 35   rhs(j) = rhs(j) - omega*t(j-1,imin)*rhs(j-1)/d(j-1)
  go to 20
!
!  far diagonals  (do vector computations).
!
 40   iwksp(imin) = iwksp(imin) + ndel
  iend = min (ibeg+ndel-1,n)
  do 45 i = ibeg,iend
 45   rhs(i) = rhs(i) - omega*t(i-ndel,imin)*rhs(i-ndel)/d(i-ndel)
  go to 20
!
!  rowwise algorithm.
!
 50   do 65 i = 1,n
     do 55 j = 1,maxt
 55      iwksp(j) = min (n,i+jt(j))
     term = omega*rhs(i)/d(i)
     do 60 j = 1,maxt
 60      rhs(iwksp(j)) = rhs(iwksp(j)) - t(i,j)*term
 65   continue
!
!  unew = inv(d)*rhs
!
 70   do 75 i = 1,n
 75   unew(i) = rhs(i)/d(i)
  return
end
subroutine sorp (ndim,nn,maxt,maxb,jt,jb,d,t,b,omega,u,rhs,unew)
!
!*******************************************************************************
!
!! SORP does an SOR solve.
!     (natural ordering, Purdue storage).
!
!        unew = inv((1/w)*d + l)*(((1-w)/w)*d*un + (rhs - u*un))
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the upper
!                triangle of the matrix
!        b      array of active size n by maxb giving the lower
!                triangle of the matrix
!        omega  over-relaxation factor
!        u      current solution vector
!        rhs    right hand side
!        unew   updated solution vector
!
!  
!
  dimension d(1), t(ndim,1), b(ndim,1), u(1), rhs(1), unew(1)
  integer   jt(ndim,1), jb(ndim,1)
  n = nn
!
!  rhs = ((1-w)/w)*d*un + (rhs - u*un)
!
  call vsubp (ndim,ndim,n,maxt,t,jt,rhs,u,unew)
  con = (1.0 - omega)/omega
  do 10 i = 1,n
 10   unew(i) = con*d(i)*u(i) + rhs(i)
!
!  unew = inv((1/w)*d + l)*rhs
!
  if (maxb >= 1) go to 20
  do 15 i = 1,n
 15   unew(i) = omega*unew(i)/d(i)
  return
 20   do 30 i = 1,n
     sum = unew(i)
     do 25 j = 1,maxb
        sum = sum - b(i,j)*unew(jb(i,j))
 25      continue
     unew(i) = omega*sum/d(i)
 30   continue
  return
end
subroutine sorstp (n,u,ubar,dnrm,ccon)
!
!*******************************************************************************
!
!! SORSTP tests if the SOR method has converged.
!
!
!  Parameters:
!
!          n      order of system
!          u      present solution estimate
!          ubar   exact solution
!          dnrm   inner product of pseudo-residuals at preceding
!                    iteration
!          con    stopping test parameter (= ccon)
!
!  
!
  dimension u(1), ubar(1)
  logical q1
  save    q1
!
! *** begin -- itpack common
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
!
! *** end   -- itpack common
!
  con = ccon
  halt = .false.
  if (ntest == 6) go to 25
!
!  special procedure for zeroth iteration.
!
  if (in >= 1) go to 5
  q1 = .false.
  udnm = 1.0
  stptst = 1000.0
  return
!
!  test if udnm needs to be recomputed
!
 5    if (q1) go to 15
  if ((in > 5)  .and.  (mod(in,5) /= 0)) go to 15
  uold = udnm
  udnm = 0.0
  do 10 i = 1,n
 10   udnm = udnm + u(i)*u(i)
  if (udnm == 0.0) udnm = 1.0
  if ((in > 5) .and.(abs (udnm-uold) <= udnm*zeta)) q1 = .true.
!
!  compute stopping test
!
 15   tr = sqrt (udnm)
  tl = 1.0
  if (con == 1.0)  go to 20
  tl = sqrt (dnrm)
  tr = tr*(1.0 - con)
 20   stptst = tl/tr
  if (tl >= tr*zeta) return
  halt = .true.
  return
!
!  second test.
!
 25   if (in == 0) ubarnm = sqrt (vdot(n,ubar,ubar))
  sum = 0.0
  do 30 i = 1,n
 30   sum = sum + (u(i) - ubar(i))**2
  tl = sqrt (sum)
  tr = ubarnm
  stptst = tl/tr
  if (tl < tr*zeta) halt = .true.
  return
end
subroutine sorw (suba,subq,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs,wksp,nw, &
  iparm,rparm,ier)
!
!*******************************************************************************
!
!! SORW drives the successive over-relaxation algorithm.
!
!
!  Parameters:
!
!          suba   matrix-vector multiplication routine
!          subq   routine to do an sor pass
!          n      input integer.  order of the system (= nn)
!          u      input/output vector.  on input, u contains the
!                 initial guess to the solution.  on output, it
!                 contains the latest estimate to the solution.
!          ubar   input vector containing the true solution
!                  (optional)
!          rhs    input vector.  contains the right hand side
!                 of the matrix problem.
!          wksp   vector used for working space.
!          nw     length of wksp array.  if this length is less than
!                  the amount needed, nw will give the needed amount
!                  upon output.
!          iparm  integer vector of length 30.  allows user to
!                 specify some integer parameters which affect
!                 the method.
!          rparm  real vector of length 30.  allows user to
!                 specify some real parameters which affect
!                 the method.
!          ier    output integer.  error flag.
!
!  
!
  external  suba, subq
  integer   iparm(30), jcoef(2), jwfac(1)
  dimension rhs(1), u(1), ubar(1), wksp(1), rparm(30), coef(1), wfac(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
!  initialize common blocks
!
  ier = 0
  n = nn
  t1 = timer (dummy)
  iacel = 3
  timit = 0.0
  digit1 = 0.0
  digit2 = 0.0
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 35
  if (level >= 2) write (nout,10)
 10   format (1x,'sor')
!
!  compute workspace base addresses and check for sufficient
!  workspace.
!
  nwksp = 2*n
  if (nw >= nwksp) go to 15
  ier = -2
  call ershow (ier,'sorw')
  go to 30
!
!  zero out workspace
!
 15   call vfill (nwksp,wksp,0.0)
!
!  iteration sequence
!
  call itsor (subq,coef,jcoef,wfac,jwfac,n,u,ubar,rhs,wksp,ier)
!
  if (ier < 0  .or.  ier == 1) go to 25
!
!  method has converged
!
  if (level >= 1) write (nout,20) in
 20   format (/1x,'sor  has converged in ',i5,' iterations' )
!
!  optional error analysis
!
 25   if (idgts < 0) go to 30
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wksp,digit1,digit2,idgts)
!
!  set return parameters in iparm and rparm
!
 30   t2 = timer (dummy)
  timit = t2 - t1
  nw = 2*n
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
  rparm(9) = omega
  rparm(10) = alphab
  rparm(11) = betab
  rparm(12) = specr
!
 35   continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
!
  return
end
subroutine split (accel,suba,subat,subq,subqt,subql,subqlt,subqr,subqrt, &
  subadp,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SPLIT determines how to apply the splitting based on IQLR.
!
  external accel, suba, subat, subq, subqt, subql, subqlt
  external subqr, subqrt, subadp, copy
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
!
!
!
  if (iqlr == 0) then
     call accel (suba,subat,copy,copy,copy,copy,subadp,coef,jcoef,n,u,ubar, &
       rhs,wksp,iwksp,iparm,rparm,jer)
  end if

  if (iqlr == 1) then
     call accel (suba,subat,subq,subqt,copy,copy,subadp,coef,jcoef,n,u, &
       ubar,rhs,wksp,iwksp,iparm,rparm,jer)
  end if

  if (iqlr == 2) then
     call accel (suba,subat,copy,copy,subq,subqt,subadp,coef,jcoef,n,u,ubar, &
       rhs,wksp,iwksp,iparm,rparm,jer)
  end if

  if (iqlr == 3) then
     call accel (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
       ubar,rhs,wksp,iwksp,iparm,rparm,jer)
  end if

  if (jer /= 0) ier = jer

  return
end
subroutine srbs (ndim,nn,maxtt,jt,d,t,omega,irwise,iwksp,x)
!
!*******************************************************************************
!
!! SRBS does an SOR back solve (natural ordering, diagonal storage).
!
!        (i + omega*inv(d)*t)*x = y
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        x      on input, x contains y
!               on output, x is the solution to back-solve
!
!  
!
  dimension x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  if (maxt <= 0) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 60
!
!  diagonal-wise algorithm.
!
  do 20 i = 1,maxt
 20   iwksp(i) = n - jt(i)
!
!  determine nc, imax.
!
 25   nc = 1
  do 30 i = 1,maxt
     nterm = iwksp(i) + 1
     if (nterm <= nc) go to 30
     nc = nterm
     imax = i
 30   continue
  if (nc <= 1) return
  ndel = jt(imax)
  iend = nc - 1
  if (ndel > 1) go to 50
!
!  special case for first super diagonal.
!
  nc1 = 1
  do 40 i = 1,maxt
     if (i == imax) go to 40
     if (iwksp(i) > nc1) nc1 = iwksp(i)
 40   continue
  iwksp(imax) = nc1 - 1
  do 45 k = iend,nc1,-1
 45   x(k) = x(k) - omega*t(k,imax)*x(k+1)/d(k)
  go to 25
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imax) = iwksp(imax) - ndel
  ibeg = max (iend - ndel,0) + 1
  do 55 i = ibeg,iend
 55   x(i) = x(i) - omega*t(i,imax)*x(i+ndel)/d(i)
  go to 25
!
!  rowwise algorithm.
!
 60   do 75 i = n,1,-1
     do 65 j = 1,maxt
 65      iwksp(j) = min (n,i+jt(j))
     sum = 0.0
     do 70 j = 1,maxt
 70      sum = sum + t(i,j)*x(iwksp(j))
     x(i) = x(i) - omega*sum/d(i)
 75   continue
  return
end
subroutine srbscp (ndim,n,jc,d,c,ncolor,nc,nt,omega,wksp,x)
!
!*******************************************************************************
!
!! SRBSCP does a back SOR solve. (Purdue storage, multicolor)
!
!     ((1/w)*d + t)*x = y
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!          x      on input, x contains y
!                 on output, x is the solution to back-solve
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1)
  dimension d(1), c(ndim,1), x(1), wksp(1)
!
  ied = n
  do 20 icol = ncolor,1,-1
     npt = nc(icol)
     ist = ied - npt + 1
     j2 = nt(icol)
     call vsubp (ndim,ndim,npt,j2,c(ist,1),jc(ist,1),x(ist),x,wksp)
     do 15 i = ist,ied
 15      x(i) = omega*x(i)/d(i)
     ied = ied - npt
 20   continue
  return
end
subroutine srbsct (ndim,n,jc,d,c,ncolor,nc,nt,nb,omega,wksp,x)
!
!*******************************************************************************
!
!! SRBSCT does a transpose back SOR solve.
!     (Purdue storage, multicolor)
!
!     ((1/w)*d + (b**t))*x = y
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length max(nc(i))
!          x      on input, x contains y
!                 on output, x is the solution to back-solve
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndim,1), x(1), wksp(1)
!
  ied = n
  do 20 icol = ncolor,1,-1
     npt = nc(icol)
     ist = ied - npt + 1
     do 15 i = ist,ied
 15      x(i) = omega*x(i)/d(i)
     j1 = nt(icol) + 1
     mj = nb(icol)
     call vsubpt (ndim,ndim,npt,mj,c(ist,j1),jc(ist,j1),x,x(ist),wksp)
     ied = ied - npt
 20   continue
  return
end
subroutine srbsp (ndim,nn,maxt,jt,d,t,omega,x)
!
!*******************************************************************************
!
!! SRBSP does an SOR backward solve (natural ordering, Purdue storage).
!        ((1/omega)*d + t)*x = y
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the upper
!                triangle of the matrix
!        omega  relaxation factor
!        x      on input, x contains y
!               on output, x is the solution to backward-solve
!
!  
!
  dimension x(1), d(1), t(ndim,1)
  integer   jt(ndim,1)
!
  n = nn
  if (maxt >= 1) go to 15
  do 10 i = 1,n
 10   x(i) = omega*x(i)/d(i)
  return
 15   do 30 i = n,1,-1
     sum = x(i)
     do 25 j = 1,maxt
        sum = sum - t(i,j)*x(jt(i,j))
 25      continue
     x(i) = omega*sum/d(i)
 30   continue
  return
end
subroutine srbst (ndim,nn,maxbb,jb,d,b,omega,irwise,iwksp,x)
!
!*******************************************************************************
!
!! SRBST does an SOR transpose back solve (natural ordering, diagonal storage).
!
!        (i + omega*inv(d)*(b**t))*x = y
!
!  Parameters:
!
!        ndim   row dimension of b array
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the sub-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxb
!        x      on input, x contains y
!               on output, x is the solution to back-solve
!
!  
!
  dimension x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
!
  n = nn
  maxb = maxbb
  if (maxb < 1) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 70
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxb
 15   iwksp(i) = n + jb(i)
!
!  determine nc, imax.
!
 20   nc = 1
  do 25 i = 1,maxb
     nterm = iwksp(i) + 1
     if (nterm <= nc) go to 25
     nc = nterm
     imax = i
 25   continue
  if (nc <= 1) return
  ndel = -jb(imax)
  iend = nc - 1
  if (ndel > 1) go to 50
!
!  special case for first sub diagonal.
!
  nc1 = 1
  do 30 i = 1,maxb
     if (i == imax) go to 30
     if (iwksp(i) > nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imax) = nc1 - 1
  do 45 k = iend,nc1,-1
 45   x(k) = x(k) - omega*b(k+1,imax)*x(k+1)/d(k)
  go to 20
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imax) = iwksp(imax) - ndel
  ibeg = max (iend - ndel,0) + 1
  do 65 i = ibeg,iend
 65   x(i) = x(i) - omega*b(i+ndel,imax)*x(i+ndel)/d(i)
  go to 20
!
!  rowwise algorithm.
!
 70   do 85 i = n,2,-1
     do 75 j = 1,maxb
 75      iwksp(j) = max (1,i+jb(j))
     term = omega*x(i)
     do 80 j = 1,maxb
 80      x(iwksp(j)) = x(iwksp(j)) - b(i,j)*term/d(iwksp(j))
 85   continue
  return
end
subroutine srbstp (ndim,nn,maxb,jb,d,b,omega,x)
!
!*******************************************************************************
!
!! SRBSTP does an SOR transpose back solve (natural ordering, Purdue storage).
!
!        ((1/omega)*d + (b**t))*x = y
!
!  Parameters:
!
!        ndim   row dimension of b array
!        n      order of system
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the lower
!                triangle of the matrix
!        omega  over-relaxation factor
!        x      on input, x contains y
!        x      on output, x is the solution to back-solve
!
!  
!
  dimension x(1), d(1), b(ndim,1)
  integer   jb(ndim,1)
!
  n = nn
  if (maxb >= 1) go to 15
  do 10 i = 1,n
 10   x(i) = omega*x(i)/d(i)
  return
 15   do 30 i = n,1,-1
     x(i) = omega*x(i)/d(i)
     term = x(i)
     do 25 j = 1,maxb
        x(jb(i,j)) = x(jb(i,j)) - b(i,j)*term
 25      continue
 30   continue
  return
end
subroutine srcg (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SRCG is the user interface to the SSOR conjugate gradient algorithm.
!
  external suba, subat, subql, subqlt, subqr, subqrt, subadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  ier = 0
  call needw ('srcg',0,irpnt,3*n+2*itmax,ier)
  if (ier < 0) return
  nw = lenr - irpnt + 1
  call srcgw (suba,subql,subadp,coef,jcoef,wksp,iwksp, &
    n,u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = irpnt + nw - 1
  return
end
subroutine srcgw (suba,subq,subadp,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs, &
  wksp,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SRCGW drives the SSOR conjugate gradient algorithm.
!
!  Parameters:
!
!          suba   matrix-vector multiplication routine
!          subq   preconditioning routine
!          subadp adpation routine
!          n      input integer.  order of the system (= nn)
!          u      input/output vector.  on input, u contains the
!                 initial guess to the solution.  on output, it
!                 contains the latest estimate to the solution.
!          ubar   input vector containing the true solution
!                  (optional)
!          rhs    input vector.  contains the right hand side
!                 of the matrix problem.
!          wksp   vector used for working space.
!          nw     length of wksp array.  if this length is less than
!                  the amount needed, nw will give the needed amount
!                  upon output.
!          iparm  integer vector of length 30.  allows user to
!                 specify some integer parameters which affect
!                 the method.
!          rparm  real vector of length 30.  allows user to
!                 specify some real parameters which affect
!                 the method.
!          ier    output integer.  error flag.
!
!  
!
  external  suba, subq, subadp
  integer   iparm(30), jcoef(2), jwfac(1)
  dimension rhs(1), u(1), ubar(1), wksp(1), rparm(30), coef(1), wfac(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
!  initialize common blocks
!
  ier = 0
  n = nn
  t1 = timer (dummy)
  iacel = 1
  timit = 0.0
  digit1 = 0.0
  digit2 = 0.0
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 35
  if (level >= 2) write (nout,10)
 10   format (1x,'srcg')
!
!  compute workspace base addresses and check for sufficient
!  workspace.
!
  iw1 = 1
  iw2 = iw1 + n
  iw3 = iw2 + n
  iw4 = iw3 + n
  nwksp = 3*n + 2*itmax
  if (nw >= nwksp) go to 15
  ier = -2
  call ershow (ier,'srcgw')
  go to 30
 15   continue
!
!  zero out workspace
!
  call vfill (nwksp,wksp,0.0)
!
!  iteration sequence
!
  call itsrcg (suba,subq,subadp,coef,jcoef,wfac,jwfac,n,u,ubar,rhs,wksp(iw1), &
    wksp(iw2),wksp(iw3),wksp(iw4),ier)
!
  if (ier < 0  .or.  ier == 1) go to 25
!
!  method has converged
!
  if (level >= 1) write (nout,20) in
 20   format (/1x,'srcg has converged in ',i5,' iterations' )
!
!  optional error analysis
!
 25   if (idgts < 0) go to 30
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wksp,digit1,digit2,idgts)
!
!  set return parameters in iparm and rparm
!
 30   t2 = timer (dummy)
  timit = t2 - t1
  nw = 3*n + 2*in
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
  rparm(9) = omega
  rparm(10) = alphab
  rparm(11) = betab
  rparm(12) = specr
!
 35   continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
!
  return
end
subroutine srfs (ndim,nn,maxbb,jb,d,b,omega,irwise,iwksp,x)
!
!*******************************************************************************
!
!! SRFS does an SOR forward solve (natural ordering, diagonal storage).
!
!        (i + omega*b*inv(d))*x = y
!
!  Parameters:
!
!        ndim   row dimension of b array
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the sub-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxb
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
!
!
  n = nn
  maxb = maxbb
  if (maxb <= 0) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 60
!
!  diagonal-wise algorithm.
!
  do 20 i = 1,maxb
 20   iwksp(i) = 1 - jb(i)
!
!  determine nc, imin.
!
 25   nc = n
  do 30 i = 1,maxb
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 30
     nc = nterm
     imin = i
 30   continue
  if (nc >= n) return
  ndel = -jb(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 50
!
!  special case for first minor subdiagonal.
!
  nc1 = n
  do 40 i = 1,maxb
     if (i == imin) go to 40
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 40   continue
  iwksp(imin) = nc1 + 1
  do 45 j = ibeg,nc1
 45   x(j) = x(j) - omega*b(j,imin)*x(j-1)/d(j-1)
  go to 25
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imin) = iwksp(imin) + ndel
  iend = min (ibeg+ndel-1,n)
  do 55 i = ibeg,iend
 55   x(i) = x(i) - omega*b(i,imin)*x(i-ndel)/d(i-ndel)
  go to 25
!
!  rowwise algorithm.
!
 60   do 75 i = 1,n
     do 65 j = 1,maxb
 65      iwksp(j) = max (1,i+jb(j))
     sum = 0.0
     do 70 j = 1,maxb
 70      sum = sum + b(i,j)*x(iwksp(j))/d(iwksp(j))
     x(i) = x(i) - omega*sum
 75   continue
  return
end
subroutine srfscp (ndim,jc,d,c,ncolor,nc,nt,nb,omega,wksp,x)
!
!*******************************************************************************
!
!! SRFSCP does a forward SOR solve.  (Purdue storage, multicolor)
!
!     ((1/w)*d + b)*x = y
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!          x      on input, x contains y
!                 on output, x is the solution to the forward solve
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndim,1), x(1), wksp(1)
!
  ist = 1
  do 20 icol = 1,ncolor
     npt = nc(icol)
     ied = ist + npt - 1
     j1 = nt(icol) + 1
     mj = nb(icol)
     call vsubp (ndim,ndim,npt,mj,c(ist,j1),jc(ist,j1),x(ist),x,wksp)
     do 15 i = ist,ied
 15      x(i) = omega*x(i)/d(i)
     ist = ist + npt
 20   continue
  return
end
subroutine srfsct (ndim,jc,d,c,ncolor,nc,nt,omega,wksp,x)
!
!*******************************************************************************
!
!! SRFSCT does a transpose forward SOR solve.  (Purdue storage, multicolor)
!
!     ((1/w)*d + (t**t))*x = y
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length max(nc(i))
!          x      on input, x contains y
!                 on output, x is the solution to the forward solve
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1)
  dimension d(1), c(ndim,1), x(1), wksp(1)
!
  ist =  1
  do 20 icol = 1,ncolor
     npt = nc(icol)
     ied = ist + npt - 1
     do 15 i = ist,ied
 15      x(i) = omega*x(i)/d(i)
     j2 = nt(icol)
     call vsubpt (ndim,ndim,npt,j2,c(ist,1),jc(ist,1),x,x(ist),wksp)
     ist = ist + npt
 20   continue
  return
end
subroutine srfsp (ndim,nn,maxb,jb,d,b,omega,x)
!
!*******************************************************************************
!
!! SRFSP does an SOR forward solve (natural ordering, Purdue storage).
!
!        ((1/omega)*d + b)*x = y
!
!  Parameters:
!
!        ndim   row dimension of b array
!        n      order of system
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the lower
!                triangle of the matrix
!        omega  relaxation factor
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), b(ndim,1)
  integer   jb(ndim,1)
!
  n = nn
  if (maxb >= 1) go to 15
  do 10 i = 1,n
 10   x(i) = omega*x(i)/d(i)
  return
 15   do 30 i = 1,n
     sum = x(i)
     do 25 j = 1,maxb
        sum = sum - b(i,j)*x(jb(i,j))
 25      continue
     x(i) = omega*sum/d(i)
 30   continue
  return
end
subroutine srfst (ndim,nn,maxtt,jt,d,t,omega,irwise,iwksp,x)
!
!*******************************************************************************
!
!! SRFST does an SOR transpose forward solve (natural ordering, diagonal storage).
!
!        (i + omega*(t**t)*inv(d))*x = y
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  n = nn
  maxt = maxtt
  if (maxt < 1) return
!
!  select rowwise or diagonal-wise algorithm.
!
  if (irwise == 1) go to 70
!
!  diagonal-wise algorithm.
!
  do 15 i = 1,maxt
 15   iwksp(i) = jt(i) + 1
!
!  determine nc, imin.
!
 20   nc = n
  do 25 i = 1,maxt
     nterm = iwksp(i) - 1
     if (nterm >= nc) go to 25
     nc = nterm
     imin = i
 25   continue
  if (nc >= n) return
  ndel = jt(imin)
  ibeg = nc + 1
  if (ndel > 1) go to 50
!
!  special case for first minor subdiagonal.
!
  nc1 = n
  do 30 i = 1,maxt
     if (i == imin) go to 30
     if (iwksp(i) < nc1) nc1 = iwksp(i)
 30   continue
  iwksp(imin) = nc1 + 1
  do 45 j = ibeg,nc1
 45   x(j) = x(j) - omega*t(j-1,imin)*x(j-1)/d(j-1)
  go to 20
!
!  far diagonals  (do vector computations).
!
 50   iwksp(imin) = iwksp(imin) + ndel
  iend = min (ibeg+ndel-1,n)
  do 65 i = ibeg,iend
 65   x(i) = x(i) - omega*t(i-ndel,imin)*x(i-ndel)/d(i-ndel)
  go to 20
!
!  rowwise algorithm.
!
 70   do 85 i = 1,n
     do 75 j = 1,maxt
 75      iwksp(j) = min (n,i+jt(j))
     term = omega*x(i)/d(i)
     do 80 j = 1,maxt
 80      x(iwksp(j)) = x(iwksp(j)) - t(i,j)*term
 85   continue
  return
end
subroutine srfstp (ndim,n,maxt,jt,d,t,omega,x)
!
!*******************************************************************************
!
!! SRFSTP does an SOR transpose forward solve (natural ordering, Purdue storage).
!
!        ((1/omega)*d + (t**t))*x = y
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the upper
!                triangle of the matrix
!        omega  over-relaxation factor
!        x      on input, x contains y
!               on output, x is the solution to forward-solve
!
!  
!
  dimension x(1), d(1), t(ndim,1)
  integer   jt(ndim,1)
!
  if (maxt >= 1) go to 15
  do 10 i = 1,n
 10   x(i) = omega*x(i)/d(i)
  return
 15   do 30 i = 1,n
     x(i) = omega*x(i)/d(i)
     term = x(i)
     do 25 j = 1,maxt
        x(jt(i,j)) = x(jt(i,j)) - t(i,j)*term
 25      continue
 30   continue
  return
end
subroutine srs (ndim,nn,maxtt,jt,d,t,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRS does an SSOR solution (natural ordering, symmetric diagonal storage).
!
!        con*(i + w*(t**t)*inv(d))*d*(i + w*inv(d)*t)*x = y
!         con = 1/(w*(2-w))   and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  fac = omega*(2.0 - omega)
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfst (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = fac*x(i)/d(i)
  call srbs (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  return
end
subroutine srs1 (ndim,nn,maxtt,jt,d,t,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRS1 does an SSOR forward solve (natural ordering, symmetric diagonal storage).
!
!        con*(i + w*(t**t)*inv(d))*d*x = y
!         con = 1/(w*(2-w))   and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  fac = omega*(2.0 - omega)
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfst (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = fac*x(i)/d(i)
  return
end
subroutine srs2 (ndim,nn,maxtt,jt,d,t,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRS2 does an SSOR back solve (natural ordering, symmetric diagonal storage).
!
!        (i + w*inv(d)*t)*x = y
!            w = omega
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  do 10 i = 1,n
 10   x(i) = y(i)
  call srbs (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  return
end
subroutine srs3 (ndim,nn,maxtt,jt,d,t,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRS3 does an SSOR transpose forward solve (natural ordering, symmetric diagonal storage).
!
!        con*d*(i + w*inv(d)*t)*x = y
!         con = 1/(w*(2-w))   and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  fac = omega*(2.0 - omega)
  do 10 i = 1,n
 10   x(i) = fac*y(i)/d(i)
  call srbs (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  return
end
subroutine srs4 (ndim,nn,maxtt,jt,d,t,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRS4 does an SSOR transpose back solve (natural ordering, symmetric diagonal storage).
!
!        (i + w*(t**t)*inv(d))*x = y
!            w = omega
!
!  Parameters:
!
!        ndim   row dimension of t array
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfst (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  return
end
subroutine srscp (ndim,nn,jc,d,c,ncolor,nc,nt,nb,omega,wksp,y,x)
!
!*******************************************************************************
!
!! SRSCP does an SSOR solve. (Purdue storage, multicolor)
!
!        con*((1/w)*d + b)*inv(d)*((1/w)*d + t)*x = y
!        where con = w/(2-w) and w = omega
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndim,1), x(1), y(1), wksp(1)
!
!
  n = nn
  fac = (2.0 - omega)/omega
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfscp (ndim,jc,d,c,ncolor,nc,nt,nb,omega,wksp,x)
  do 15 i = 1,n
 15   x(i) = fac*d(i)*x(i)
  call srbscp (ndim,n,jc,d,c,ncolor,nc,nt,omega,wksp,x)
  return
end
subroutine srscp1 (ndim,nn,jc,d,c,ncolor,nc,nt,nb,omega,wksp,y,x)
!
!*******************************************************************************
!
!! SRSCP1 does an SSOR forward solve.  (Purdue storage, multicolor)
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndim,1), x(1), y(1), wksp(1)
!
!
  n = nn
  fac = (2.0 - omega)/omega
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfscp (ndim,jc,d,c,ncolor,nc,nt,nb,omega,wksp,x)
  do 15 i = 1,n
 15   x(i) = fac*d(i)*x(i)
  return
end
subroutine srscp2 (ndim,n,jc,d,c,ncolor,nc,nt,omega,wksp,y,x)
!
!*******************************************************************************
!
!! SRSCP2 does an SSOR back solve.  (Purdue storage, multicolor)
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1)
  dimension d(1), c(ndim,1), x(1), y(1), wksp(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srbscp (ndim,n,jc,d,c,ncolor,nc,nt,omega,wksp,x)
  return
end
subroutine srscp3 (ndim,n,jc,d,c,ncolor,nc,nt,nb,omega,wksp,y,x)
!
!*******************************************************************************
!
!! SRSCP3 does a transpose SSOR back solve.  (Purdue storage, multicolor)
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length max(nc(i))
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndim,1), x(1), y(1), wksp(1)
!
!
  fac = (2.0 - omega)/omega
!
  do 15 i = 1,n
 15   x(i) = fac*d(i)*y(i)
  call srbsct (ndim,n,jc,d,c,ncolor,nc,nt,nb,omega,wksp,x)
  return
end
subroutine srscp4 (ndim,n,jc,d,c,ncolor,nc,nt,omega,wksp,y,x)
!
!*******************************************************************************
!
!! SRSCP4 does a transpose ssor forward solve.  (Purdue storage, multicolor)
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length max(nc(i))
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1)
  dimension d(1), c(ndim,1), x(1), y(1), wksp(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfsct (ndim,jc,d,c,ncolor,nc,nt,omega,wksp,x)
  return
end
subroutine srscpt (ndim,nn,jc,d,c,ncolor,nc,nt,nb,omega,wksp,y,x)
!
!*******************************************************************************
!
!! SRSCPT does a transpose SSOR solve.  (Purdue storage, multicolor)
!
!  Parameters:
!
!          ndim   row dimension of c,jc arrays
!          n      order of system (= nn)
!          jc     integer array giving the column indices of the
!                  corresponding elements in c
!          d      vector of length n giving the diagonal elements
!                  of the matrix
!          c      array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!          ncolor number of colors used
!          nc     integer vector of length ncolor giving the number
!                  of nodes for each color
!          nt     integer vector of length ncolor giving the number
!                  of upper columns for each color
!          nb     integer vector of length ncolor giving the number
!                  of lower columns for each color
!          omega  over-relaxation factor
!          wksp   workspace vector of length max(nc(i))
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndim,1), x(1), y(1), wksp(1)
!
!
  n = nn
  fac = (2.0 - omega)/omega
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfsct (ndim,jc,d,c,ncolor,nc,nt,omega,wksp,x)
  do 15 i = 1,n
 15   x(i) = fac*d(i)*x(i)
  call srbsct (ndim,n,jc,d,c,ncolor,nc,nt,nb,omega,wksp,x)
  return
end
subroutine srsi (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n,u, &
  ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SRSI is the user interface to the SSOR Chebyshev acceleration algorithm.
!
  external suba, subat, subql, subqlt, subqr, subqrt, subadp
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  ier = 0
  call needw ('srsi',0,irpnt,4*n,ier)
  if (ier < 0) return
  nw = lenr - irpnt + 1
  call srsiw (suba,subql,subadp,coef,jcoef,wksp,iwksp, &
    n,u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = irpnt + nw - 1
  return
end
subroutine srsiw (suba,subq,subadp,coef,jcoef,wfac,jwfac,nn,u,ubar,rhs,wksp, &
  nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SRSIW drives the SSOR Chebyshev acceleration algorithm.
!
!  Parameters:
!
!          suba   matrix-vector multiplication routine
!          subq   preconditioning routine
!          subadp adpation routine
!          n      input integer.  order of the system (= nn)
!          u      input/output vector.  on input, u contains the
!                 initial guess to the solution.  on output, it
!                 contains the latest estimate to the solution.
!          ubar   input vector containing the true solution
!                  (optional)
!          rhs    input vector.  contains the right hand side
!                 of the matrix problem.
!          wksp   vector used for working space.
!          nw     length of wksp array.  if this length is less than
!                  the amount needed, nw will give the needed amount
!                  upon output.
!          iparm  integer vector of length 30.  allows user to
!                 specify some integer parameters which affect
!                 the method.
!          rparm  real vector of length 30.  allows user to
!                 specify some real parameters which affect
!                 the method.
!          ier    output integer.  error flag.
!
!  
!
  external  suba, subq, subadp
  integer   iparm(30), jcoef(2), jwfac(1)
  dimension rhs(1), u(1), ubar(1), wksp(1), rparm(30), coef(1),wfac(1)
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
!  initialize common blocks
!
  ier = 0
  n = nn
  t1 = timer (dummy)
  iacel = 2
  timit = 0.0
  digit1 = 0.0
  digit2 = 0.0
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 35
  if (level >= 2) write (nout,10)
 10   format (1x,'srsi')
!
!  compute workspace base addresses and check for sufficient
!  workspace.
!
  iw1 = 1
  iw2 = iw1 + n
  iw3 = iw2 + n
  iw4 = iw3 + n
  nwksp = 4*n
  if (nw >= nwksp) go to 15
  ier = -2
  call ershow (ier,'srsiw')
  go to 30
 15   continue
!
!  compute an initial rayleigh quotient and adjust emax, emin.
!
  call vfill (n,wksp,1.0)
  call subq (coef,jcoef,wfac,jwfac,n,wksp,wksp(iw2))
  call suba (coef,jcoef,wfac,jwfac,n,wksp(iw2),wksp(iw3))
  rq = vdot (n,wksp(iw2),wksp(iw3)) / vdot (n,wksp(iw2),wksp)
  rqmax = 1.0
  rqmin = rq
!
!  adjust emax, emin.
!
  emax = 1.0
  maxadd = .false.
  if (minadd) emin = amin1 (emin,rqmin)
  if (minadd) emin = amax1 (emin,0.0)
!
!  zero out workspace
!
  call vfill (nwksp,wksp,0.0)
!
!  iteration sequence
!
  call itsrsi (suba,subq,subadp,coef,jcoef,wfac,jwfac,n,u,ubar, &
    rhs,wksp(iw1),wksp(iw2),wksp(iw3),wksp(iw4),ier)
!
  if (ier < 0  .or.  ier == 1) go to 25
!
!  method has converged
!
  if (level >= 1) write (nout,20) in
 20   format (/1x,'srsi  has converged in ',i5,' iterations ')
!
!  optional error analysis
!
 25   if (idgts < 0) go to 30
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wksp,digit1,digit2,idgts)
!
!  set return parameters in iparm and rparm
!
 30   t2 = timer (dummy)
  timit = t2 - t1
  nw = 4*n
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
  rparm(9) = omega
  rparm(10) = alphab
  rparm(11) = betab
  rparm(12) = specr
!
 35   continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
!
  return
end
subroutine srsn (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRSN does an SSOR solution (natural ordering, nonsymmetric diagonal storage).
!
!        con*(i + w*b*inv(d))*d*(i + w*inv(d)*t)*x = y
!         where  con = 1/(w*(2-w))  and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        b      array of active size n by maxb giving the sub-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1), b(ndim,1)
  integer   jt(1), jb(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  fac = omega*(2.0 - omega)
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfs (ndim,n,maxb,jb,d,b,omega,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = fac*x(i)/d(i)
  call srbs (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  return
end
subroutine srsn1 (ndim,n,maxb,jb,d,b,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRSN1 does an SSOR forward pass (natural ordering, nonsymmetric diagonal storage).
!
!        con*(i + w*b*inv(d))*d*(i + w*inv(d)*t)*x = y
!         where  con = 1/(w*(2-w))  and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the sub-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
!
!
  fac = omega*(2.0 - omega)
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfs (ndim,n,maxb,jb,d,b,omega,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = fac*x(i)/d(i)
  return
end
subroutine srsn2 (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRSN2 does an SSOR backward pass (natural ordering, nonsymmetric diagonal storage).
!
!        con*(i + w*b*inv(d))*d*(i + w*inv(d)*t)*x = y
!         where  con = 1/(w*(2-w))  and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srbs (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  return
end
subroutine srsn3 (ndim,n,maxb,jb,d,b,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRSN3 does a transpose SSOR back pass (natural ordering, nonsymmetric diagonal storage).
!
!       con*(i + w*(t**t)*inv(d))*d*(i + w*inv(d)*(b**t))*x = y
!        con = 1/(w*(2-w))  and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxb   number of columns in b array
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the sub-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndim,1)
  integer   jb(1), iwksp(1)
!
!
  fac = omega*(2.0 - omega)
  do 15 i = 1,n
 15   x(i) = fac*y(i)/d(i)
  call srbst (ndim,n,maxb,jb,d,b,omega,irwise,iwksp,x)
  return
end
subroutine srsn4 (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRSN4 does a transpose SSOR forward pass (natural ordering, nonsymmetric diagonal storage).
!
!       con*(i + w*(t**t)*inv(d))*d*(i + w*inv(d)*(b**t))*x = y
!        con = 1/(w*(2-w))  and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(1), iwksp(1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfst (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  return
end
subroutine srsnt (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,omega,irwise,iwksp,y,x)
!
!*******************************************************************************
!
!! SRSNT does a transpose SSOR solution (natural ordering, nonsymmetric diagonal storage).
!
!       con*(i + w*(t**t)*inv(d))*d*(i + w*inv(d)*(b**t))*x = y
!        con = 1/(w*(2-w))  and  w = omega
!
!  Parameters:
!
!        ndim   row dimension of t and b arrays
!        n      order of system (= nn)
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer vector of length maxt giving the diagonal
!                indices of the corresponding columns in t
!        jb     integer vector of length maxb giving the diagonal
!                indices of the corresponding columns in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the super-
!                diagonals of the matrix
!        b      array of active size n by maxb giving the sub-
!                diagonals of the matrix
!        omega  over-relaxation factor
!        irwise rowwise algorithm switch
!                = 0  use diagonal algorithm
!                = 1  use row-wise algorithm
!        iwksp  integer workspace of length maxt
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1), b(ndim,1)
  integer   jt(1), jb(1), iwksp(1)
!
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  fac = omega*(2.0 - omega)
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfst (ndim,n,maxt,jt,d,t,omega,irwise,iwksp,x)
  do 15 i = 1,n
 15   x(i) = fac*x(i)/d(i)
  call srbst (ndim,n,maxb,jb,d,b,omega,irwise,iwksp,x)
  return
end
subroutine srsntp (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,omega,y,x)
!
!*******************************************************************************
!
!! SRSNTP does an SSOR transpose solution (natural ordering, Purdue storage).
!
!        con*((1/w)*d + (t**t))*inv(d)*((1/w)*d + (b**t))*x = y
!        where con = w/(2-w) and w = omega
!
!  Parameters:
!
!        ndim   row dimension of t,b arrays
!        n      order of system
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the upper
!                triangle of the matrix
!        b      array of active size n by maxb giving the lower
!                triangle of the matrix
!        omega  relaxation factor
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1), b(ndim,1)
  integer   jt(ndim,1), jb(ndim,1)
!
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  fac = (2.0 - omega)/omega
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfstp (ndim,n,maxt,jt,d,t,omega,x)
  do 15 i = 1,n
 15   x(i) = fac*d(i)*x(i)
  call srbstp (ndim,n,maxb,jb,d,b,omega,x)
  return
end
subroutine srsp (ndim,nn,maxtt,maxbb,jt,jb,d,t,b,omega,y,x)
!
!*******************************************************************************
!
!! SRSP does an SSOR solution (natural ordering Purdue storage).
!
!        con*((1/w)*d + b)*inv(d)*((1/w)*d + t)*x = y
!        where con = w/(2-w) and w = omega
!
!  Parameters:
!
!        ndim   row dimension of t,b arrays
!        n      order of system
!        maxt   number of columns in t array
!        maxb   number of columns in b array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the upper
!                triangle of the matrix
!        b      array of active size n by maxb giving the lower
!                triangle of the matrix
!        omega  relaxation factor
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1), b(ndim,1)
  integer   jt(ndim,1), jb(ndim,1)
!
!
  n = nn
  maxt = maxtt
  maxb = maxbb
  fac = (2.0 - omega)/omega
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfsp (ndim,n,maxb,jb,d,b,omega,x)
  do 15 i = 1,n
 15   x(i) = fac*d(i)*x(i)
  call srbsp (ndim,n,maxt,jt,d,t,omega,x)
  return
end
subroutine srsp1 (ndim,n,maxb,jb,d,b,omega,y,x)
!
!*******************************************************************************
!
!! SRSP1 does an SSOR forward solve (natural ordering, Purdue storage).
!
!  Parameters:
!
!        ndim   row dimension of t,b arrays
!        n      order of system
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the lower
!                triangle of the matrix
!        omega  relaxation factor
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndim,1)
  integer   jb(ndim,1)
!
!
  fac = (2.0 - omega)/omega
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfsp (ndim,n,maxb,jb,d,b,omega,x)
  do 15 i = 1,n
 15   x(i) = fac*d(i)*x(i)
  return
end
subroutine srsp2 (ndim,n,maxt,jt,d,t,omega,y,x)
!
!*******************************************************************************
!
!! SRSP2 does an SSOR back solve (natural ordering, Purdue storage).
!
!  Parameters:
!
!        ndim   row dimension of t,b arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the upper
!                triangle of the matrix
!        omega  relaxation factor
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(ndim,1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srbsp (ndim,n,maxt,jt,d,t,omega,x)
  return
end
subroutine srsp3 (ndim,n,maxb,jb,d,b,omega,y,x)
!
!*******************************************************************************
!
!! SRSP3 does an SSOR transpose back solve (natural ordering, Purdue storage).
!
!  Parameters:
!
!        ndim   row dimension of t,b arrays
!        n      order of system
!        maxb   number of columns in b array
!        jb     integer array giving the column numbers of the
!                corresponding elements in b
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        b      array of active size n by maxb giving the lower
!                triangle of the matrix
!        omega  relaxation factor
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), b(ndim,1)
  integer   jb(ndim,1)
!
!
  fac = (2.0 - omega)/omega
  do 15 i = 1,n
 15   x(i) = fac*d(i)*y(i)
  call srbstp (ndim,n,maxb,jb,d,b,omega,x)
  return
end
subroutine srsp4 (ndim,n,maxt,jt,d,t,omega,y,x)
!
!*******************************************************************************
!
!! SRSP4 does an SSOR transpose forward solve (natural ordering, Purdue storage).
!
!  Parameters:
!
!        ndim   row dimension of t,b arrays
!        n      order of system
!        maxt   number of columns in t array
!        jt     integer array giving the column numbers of the
!                corresponding elements in t
!        d      vector of length n giving the diagonal elements
!                of the matrix
!        t      array of active size n by maxt giving the upper
!                triangle of the matrix
!        omega  relaxation factor
!        y      right-hand-side vector
!        x      on output, x is the solution
!
!  
!
  dimension y(1), x(1), d(1), t(ndim,1)
  integer   jt(ndim,1)
!
  do 10 i = 1,n
 10   x(i) = y(i)
  call srfstp (ndim,n,maxt,jt,d,t,omega,x)
  return
end
subroutine ssor1 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SSOR1 drives the point SSOR method.
!
  external accel, suba8, suba9, subq79, subq80, subq81, subq82
  external subq83, subq84, subq85
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  n = nn
  iwkpt1 = irpnt
  irpnt = irpnt + n
  if (isymm /= 0) irpnt = irpnt + n
  call move1 (ndim,mdim,n,maxnz,jcoef,coef,maxt,maxb,ier)
  if (ier < 0) then
     call ershow (ier,'ssor1')
     return
  end if
  call split (accel,suba8,suba9,subq79,subq80,subq81,subq82,subq83,subq84, &
    subq85,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  if (isymm /= 0) irpnt = irpnt - n
  irpnt = irpnt - n
  return
end
subroutine ssor2 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SSOR2 drives the point SSOR method.
!
  external accel, suba1, subq7, subq8, subq9, subq10,subq11, subq12
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  call rowise (maxnz,jcoef,irwise)
  call needw ('ssor2',1,iipnt,maxnz,ier)
  if (ier < 0) return
  iwkpt1 = iipnt
  iipnt = iipnt + maxnz
  call split (accel,suba1,suba1,subq7,subq7,subq8,subq9,subq10,subq11,subq12, &
    coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  iipnt = iipnt - maxnz
  return
end
subroutine ssor3 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SSOR3 drives the point SSOR method.
!
  external accel, suba4, suba5, subq41, subq42, subq43, subq44
  external subq45, subq46, subq47
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  call rowise (maxnz,jcoef,irwise)
  call needw ('ssor3',1,iipnt,maxnz,ier)
  if (ier < 0) return
  call needw ('ssor3',0,irpnt,n,ier)
  if (ier < 0) return
  call move2 (ndim,n,maxnz,jcoef,coef,wksp(irpnt),iwksp(iipnt),maxt,maxb)
  iwkpt1 = irpnt
  irpnt = irpnt + n
  iwkpt2 = iipnt
  iipnt = iipnt + maxnz
  call split (accel,suba4,suba5,subq41,subq42,subq43,subq44,subq45,subq46, &
    subq47,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  iipnt = iipnt - maxnz
  return
end
subroutine ssor6 (accel,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SSOR6 drives the multi-color SSOR method.
!
  external accel, suba8, suba9, subq97, subq98, subq99, sub100
  external sub101, sub102, sub103
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
!
  iwkpt1 = irpnt
  irpnt = irpnt + n + ncmax
  call split (accel,suba8,suba9,subq97,subq98,subq99,sub100,sub101,sub102, &
    sub103,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n - ncmax
  return
end
subroutine ssor7 (accel,coef,jcoef,nn,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! SSOR7 drives the multi-color SSOR method.
!
  external accel, suba2, suba3, subq27, subq28, subq29
  external subq30, subq31, subq32, subq33
  integer   iparm(30), jcoef(2), iwksp(1)
  dimension rhs(1), u(1), ubar(1), rparm(30), coef(1), wksp(1)
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
!
  n = nn
  t1 = timer (dummy)
  if (ifact == 1) call mfact (coef,jcoef,wksp,iwksp,n,ier)
  t2 = timer (dummy)
  timfac = t2 - t1
  if (ier < 0) return
  iwkpt1 = irpnt
  irpnt = irpnt + n
  call split (accel,suba2,suba3,subq27,subq28,subq29,subq30,subq31,subq32, &
    subq33,coef,jcoef,n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
  irpnt = irpnt - n
  return
end
subroutine ssorad (ssorcp,coef,jcoef,wfac,jwfac,n,p,z,r,icode)
!
!*******************************************************************************
!
!! SSORAD does the SSOR adaptive process.
!
!  Parameters:
!
!         n       order of system
!         p,z,r   vectors from acceleration algorithm
!         icode   key for restarting iteration
!                  = 0    omega unchanged (no restart)
!                  = 1    new omega       (restart needed)
!
!  
!
  dimension p(1), z(1), r(1), coef(1), jcoef(2), wfac(1), jwfac(1)
  external ssorcp
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
!
!
!
!  parameter estimation formulas
!
  alp (w,beta,s) = ((1.0 + beta*w*w)*s - w*(2.0 - w)) / (w*(2.0 - w - s))
!
  omg (alpha,beta) = 2.0/(1.0 + sqrt (1.0 + 2.0*alpha + 4.0*beta))
!
  se (w,alpha,beta) = ((1.0 + alpha)*w*(2.0 - w)) / (1.0 + alpha*w + beta*w*w)
!
  cond (w,alpha,beta) = 1.0/se(w,alpha,beta)
!
  rc (w,alpha,beta) = alog ((sqrt (cond(w,alpha,beta))+1.0) / &
    (sqrt (cond(w,alpha,beta))-1.0))
!
  icode = 0
  if (is >= 6  .and.  (.not. minadp)) go to 5
  tmo = 2.0 - omega
  if (emin < tmo) alphab = amin1 (alphab, alp(omega,betab,emin))
 5    if ((.not. omgadp) .or. (.not. minadp) .or. (is <= 5)) return
  omegab = amax1 (1.0, omg (alphab,betab))
  if (rc(omega,alphab,betab) > fff*rc(omegab,alphab,betab)) return
  if (iacel == 2) pap = vdot (n,p,z)
  call omgchg (ssorcp,coef,jcoef,wfac,jwfac,n,p,r)
  omega = amax1 (1.0,omg(alphab,betab))
  icode = 1
  if (level >= 2) write (nout,10) in, alphab, betab, omega
 10   format (/1x,15x,36hparameters were changed at iteration,i7/ &
     1x,20x,19halphab             ,f15.9/ &
     1x,20x,19hbetab              ,f15.9/ &
     1x,20x,19homega              ,f15.9/)
  return
end
subroutine ssord (ndim,maxt,jt,d,t,nn,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SSORD computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for symmetric diagonal storage format.
!
!  Parameters:
!
!         ndim    row dimension of coef array in defining routine
!         maxt    number of diagonals in t
!         jt      diagonal numbers for upper triangular part
!         d       diagonal
!         t       upper triangular diagonals
!         n       order of system
!         p       vector from acceleration algorithm
!         r       workspace vector from acceleration algorithm
!         pdp     (p,d*p)
!         pldup   (p,l*d*u*p)
!
!  
!
  integer   jt(1)
  dimension d(1), t(ndim,1), p(1), r(1)
!
!  compute pdp = (p,d*p).
!
  n = nn
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*d(i)*p(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p) = (u*p,inv(d)*u*p)
!
  pldup = 0.0
  if (maxt <= 0) return
  do 15 i = 1,n
 15   r(i) = 0.0
  call vaddd (ndim,1,n,n,maxt,t,jt,r,p,0)
  sum = 0.0
  do 20 i = 1,n
 20   sum = sum + r(i)*r(i)/d(i)
  pldup = sum
  return
end
subroutine ssordn (ndim,maxt,maxb,jt,jb,d,t,b,nn,p,r,wksp,pdp,pldup)
!
!*******************************************************************************
!
!! SSORDN computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for nonsymmetric diagonal storage format.
!
!  Parameters:
!
!         ndim    row dimension of coef array in defining routine
!         maxt    number of diagonals in t
!         maxb    number of diagonals in b
!         jt      diagonal numbers for upper triangular part
!         jb      diagonal numbers for lower triangular part
!         d       diagonal
!         t       upper triangular diagonals
!         b       lower triangular diagonals
!         n       order of system
!         p       vector from acceleration algorithm
!         r       workspace vector from acceleration algorithm
!         wksp    workspace vector of length n
!         pdp     (p,d*p)
!         pldup   (p,l*d*u*p)
!
!  
!
  integer   jt(1), jb(1)
  dimension d(1), t(ndim,1), b(ndim,1), p(1), r(1), wksp(1)
!
!  compute pdp = (p,d*p).
!
  n = nn
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*d(i)*p(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p)
!
  pldup = 0.0
  if (maxt <= 0 .or. maxb <= 0) return
  do 15 i = 1,n
 15   r(i) = 0.0
  call vaddd (ndim,1,n,n,maxt,t,jt,r,p,0)
  do 20 i = 1,n
 20   r(i) = r(i)/d(i)
  do 25 i = 1,n
 25   wksp(i) = 0.0
  call vaddd (ndim,1,n,n,maxb,b,jb,wksp,r,0)
  sum = 0.0
  do 30 i = 1,n
 30   sum = sum + p(i)*wksp(i)
  pldup = sum
  return
end
subroutine ssorp (ndim,maxt,jt,d,t,nn,p,r,wksp,pdp,pldup)
!
!*******************************************************************************
!
!! SSORP computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for symmetric Purdue storage format.
!
!  Parameters:
!
!         ndim    row dimension of coef array in defining routine
!         maxt    number of columns in t
!         jt      column numbers for upper triangular part
!         d       diagonal
!         t       upper triangular part of a
!         n       order of system
!         p       vector from acceleration algorithm
!         r       workspace vector from acceleration algorithm
!         wksp    workspace vector of length n
!                  (keygs = 1 only)
!         pdp     (p,d*p)
!         pldup   (p,l*d*u*p)
!
!  
!
  integer   jt(ndim,1)
  dimension d(1), t(ndim,1), p(1), r(1), wksp(1)
!
!  compute pdp = (p,d*p).
!
  n = nn
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*d(i)*p(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p) = (u*p,inv(d)*u*p)
!
  pldup = 0.0
  if (maxt <= 0) return
  do 15 i = 1,n
 15   r(i) = 0.0
  call vaddp (ndim,ndim,n,maxt,t,jt,r,p,wksp)
  sum = 0.0
  do 20 i = 1,n
 20   sum = sum + r(i)*r(i)/d(i)
  pldup = sum
  return
end
subroutine ssorpn (ndimm,maxt,maxb,jt,jb,d,t,b,nn,p,r,wksp,pdp,pldup)
!
!*******************************************************************************
!
!! SSORPN computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for nonsymmetric Purdue storage format.
!
!  Parameters:
!
!         ndim    row dimension of coef array in defining routine
!         maxt    number of columns in t
!         maxb    number of columns in b
!         jt      column numbers for upper triangular part
!         jb      column numbers for lower triangular part
!         d       diagonal
!         t       upper triangular part
!         b       lower triangular part
!         n       order of system
!         p       vector from acceleration algorithm
!         r       workspace vector from acceleration algorithm
!         wksp    workspace vector of length n
!                  2*n if keygs = 1
!         pdp     (p,d*p)
!         pldup   (p,l*d*u*p)
!
!  
!
  integer   jt(ndimm,1), jb(ndimm,1)
  dimension d(1), t(ndimm,1), b(ndimm,1), p(1), r(1), wksp(1)
!
!  compute pdp = (p,d*p).
!
  n = nn
  ndim = ndimm
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*d(i)*p(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p)
!
  pldup = 0.0
  if (maxt <= 0 .or. maxb <= 0) return
  do 15 i = 1,n
 15   r(i) = 0.0
  call vaddp (ndim,ndim,n,maxt,t,jt,r,p,wksp)
  do 20 i = 1,n
 20   r(i) = r(i)/d(i)
  do 25 i = 1,n
 25   wksp(i) = 0.0
  np1 = n + 1
  call vaddp (ndim,ndim,n,maxb,b,jb,wksp,r,wksp(np1))
  sum = 0.0
  do 30 i = 1,n
 30   sum = sum + p(i)*wksp(i)
  pldup = sum
  return
end
subroutine ssrcd (ldf,ndim,maxnz,nsize,iblock,dfac,coef,jcoef,nn,p,r, &
  wksp,pdp,pldup)
!
!*******************************************************************************
!
!! SSRCD computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for symmetric block diagonal storage format.
!
!  Parameters:
!
!         ldf      row dimension of dfac
!         ndim     row dimension of coef array
!         maxnz    number of diagonals stored in coef
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         dfac     array for diagonal block factorization
!         coef     array for matrix coefficients
!         jcoef    vector for diagonal numbers
!         n        size of system
!         p        vector from acceleration algorithm
!         r        workspace vector from acceleration algorithm
!         wksp     workspace vector of length n
!         pdp      (p,d*p)
!         pldup    (p,l*d*u*p)
!
!  
!
  integer   jcoef(2), iblock(3,1)
  dimension dfac(ldf,1), coef(ndim,2), p(1), r(1), wksp(1)
!
!  compute pdp = (p,d*p).
!
  n = nn
  nwdiag = iblock (3,1)
  nt = nwdiag - 1
  call bmul (ndim,n,nt,coef,coef(1,2),p,r)
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*r(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p) = (u*p,inv(d)*u*p)
!
  do 15 i = 1,n
 15   r(i) = 0.0
  jbgn = nwdiag + 1
  mdiag = maxnz - nwdiag
  call vaddd (ndim,1,n,n,mdiag,coef(1,jbgn),jcoef(jbgn),r,p,0)
  call bdsol (ldf,n,nsize,nt,0,dfac,r,wksp,0)
  sum = 0.0
  do 25 i = 1,n
 25   sum = sum + r(i)*wksp(i)
  pldup = sum
  return
end
subroutine ssrcdm (ldf,ndim,lbhb,nsize,ncol,nci,ipt,iblock,dfac,coef,jcnew, &
  nn,p,r,wksp,pdp,pldup)
!
!*******************************************************************************
!
!! SSRCDM computes pdp = (p,d*p)  and pldup = (p,l*inv(d)*u*p).
!
!     for nonsymmetric block diagonal storage format.
!     (nonconstant block size)
!
!  Parameters:
!
!         ldf      row dimension of dfac array
!         ndim     row dimension of coef array
!         lbhb     integer vector of size ncolor giving the number
!                   of diagonal blocks for each distinct block size.
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         ncolor   number of distinct block sizes
!         nci      integer vector of length ncolor, giving the number
!                   of nodes for each distinct block size.
!         ipt      integer pointer vector of length ncolor+1
!                   giving the starting locations of new block
!                   rows
!         iblock   integer array of size 3 by ncolor by max(lbhb(i))
!                   giving block constants
!         dfac     array for diagonal block factorization
!         coef     array of matrix coefficients
!         jcnew    integer array of row dimension ncolor giving the
!                   diagonal numbers for each block
!         n        size of system
!         p        vector from acceleration algorithm
!         r        workspace vector from acceleration algorithm
!         wksp     workspace vector of length n
!         pdp      (p,d*p)
!         pldup    (p,l*d*u*p)
!
!  
!
  integer   jcnew(ncol,1), iblock(3,ncol,2), lbhb(1), nci(1), ipt(1)
  dimension dfac(ldf,1), coef(ndim,2), p(1), r(1), wksp(1)
!
!  define constants ndt, ndb.
!
  n = nn
  ncolor = ncol
  ndt = iblock(3,1,1) - 1
  ndb = iblock(3,1,2)
  nwdiag = ndt + ndb + 1
!
!  compute pdp = (p,d*p).
!
  ind = ndt + 2
  call bmuln (ndim,n,ndt,ndb,coef,coef(1,2),coef(1,ind),p,r)
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*r(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p)
!
  do 15 i = 1,n
     r(i) = 0.0
     wksp(i) = 0.0
 15   continue
  do 25 k = 1,ncolor
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     do 20 j = 3,jlim
        jcol = k + iblock(1,k,j)
        if (jcol <= k) go to 20
        jstb = iblock(2,k,j) + nwdiag
        mb = iblock(3,k,j)
        inc = ipt(jcol) - ipt(k)
        nb = nci(jcol)
        istb = ist + inc
        call vaddd (ndim,ncolor,na,nb,mb,coef(ist,jstb),jcnew(k,jstb), &
          r(ist),p(istb),inc)
 20      continue
 25   continue
  call bdsol (ldf,n,nsize,ndt,ndb,dfac,r,r,1)
  do 35 k = 1,ncolor
     ist = ipt(k) + 1
     jlim = lbhb(k)
     na = nci(k)
     do 30 j = 3,jlim
        jcol = k + iblock(1,k,j)
        if (jcol >= k) go to 30
        jstb = iblock(2,k,j) + nwdiag
        mb = iblock(3,k,j)
        inc = ipt(jcol) - ipt(k)
        nb = nci(jcol)
        istb = ist + inc
        call vaddd (ndim,ncolor,na,nb,mb,coef(ist,jstb),jcnew(k,jstb), &
          wksp(ist),r(istb),inc)
 30      continue
 35   continue
  sum = 0.0
  do 40 i = 1,n
 40   sum = sum + p(i)*wksp(i)
  pldup = sum
  return
end
subroutine ssrcdn (ldf,ndim,lbhb,nsize,iblock,dfac,coef,jcoef,nn,p,r,wksp, &
  pdp,pldup)
!
!*******************************************************************************
!
!! SSRCDN computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for nonsymmetric block diagonal storage format.
!     (constant block size)
!
!  Parameters:
!
!         ldf      row dimension of dfac
!         ndim     row dimension of coef array
!         lbhb     number of blocks per block row
!         nsize    size of an individual subsystem within a
!                   diagonal block
!         iblock   integer array of size 3 by lbhb
!                   giving block constants
!         dfac     array for diagonal block factorization
!         coef     array for matrix coefficients
!         jcoef    vector for diagonal numbers
!         n        size of system
!         p        vector from acceleration algorithm
!         r        workspace vector from acceleration algorithm
!         wksp     workspace vector of length n
!         pdp      (p,d*p)
!         pldup    (p,l*d*u*p)
!
!  
!
  integer   jcoef(2), iblock(3,2)
  dimension dfac(ldf,1), coef(ndim,2), p(1), r(1), wksp(1)
!
!  compute nt, nb, maxt, maxb
!
  n = nn
  nt = iblock(3,1) - 1
  nb = iblock(3,2)
  maxt = 0
  maxb = 0
  if (lbhb < 3) go to 15
  do 10 j = 3,lbhb
     ind = iblock(1,j)
     if (ind > 0) maxt = maxt + iblock(3,j)
     if (ind < 0) maxb = maxb + iblock(3,j)
 10   continue
!
!  compute pdp = (p,d*p).
!
 15   ind = nt + 2
  call bmuln (ndim,n,nt,nb,coef,coef(1,2),coef(1,ind),p,r)
  sum = 0.0
  do 20 i = 1,n
 20   sum = sum + p(i)*r(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p)
!
  do 25 i = 1,n
     wksp(i) = 0.0
     r(i) = 0.0
 25   continue
  ind = nt + nb + 2
  indd = ind + maxt
  call vaddd (ndim,1,n,n,maxt,coef(1,ind),jcoef(ind),r,p,0)
  call bdsol (ldf,n,nsize,nt,nb,dfac,r,r,1)
  call vaddd (ndim,1,n,n,maxb,coef(1,indd),jcoef(indd),wksp,r,0)
  sum = 0.0
  do 30 i = 1,n
 30   sum = sum + p(i)*wksp(i)
  pldup = sum
  return
end
subroutine ssrcp (ndim,jc,d,c,nn,ncolor,nc,nt,p,r,wksp,pdp,pldup)
!
!*******************************************************************************
!
!! SSRCP computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for symmetric multicolor Purdue storage format.
!
!  Parameters:
!
!         ndim    row dimension of c,jc arrays
!         jc      integer array giving the column indices of the
!                  corresponding elements in c
!         d       vector of length n giving the diagonal elements
!                  of the matrix
!         c       array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!         n       order of system
!         ncolor  number of colors used
!         nc      integer vector of length ncolor giving the number
!                  of nodes for each color
!         nt      integer vector of length ncolor giving the number
!                  of upper columns for each color
!         p       vector from acceleration algorithm
!         r       workspace vector from acceleration algorithm
!         wksp    workspace vector of length
!                  max(nc(i))     if keygs = 1
!                  0              if keygs = 2
!         pdp     (p,d*p)
!         pldup   (p,l*d*u*p)
!
!  
!
  integer   jc(ndim,1), nc(1), nt(1)
  dimension d(1), c(ndim,1), p(1), r(1), wksp(1)
!
!  compute pdp = (p,d*p).
!
  n = nn
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*d(i)*p(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p) = (u*p,inv(d)*u*p)
!
  do 15 i = 1,n
 15   r(i) = 0.0
  ist = 1
  do 20 icol = 1,ncolor
     npt = nc(icol)
     mj = nt(icol)
     call vaddp (ndim,ndim,npt,mj,c(ist,1),jc(ist,1),r(ist),p,wksp)
     ist = ist + npt
 20   continue
  sum = 0.0
  do 25 i = 1,n
 25   sum = sum + r(i)*r(i)/d(i)
  pldup = sum
  return
end
subroutine ssrcpn (ndimm,jc,d,c,nn,ncol,nc,nt,nb,p,r,wksp,pdp,pldup)
!
!*******************************************************************************
!
!! SSRCPN computes pdp = (p,d*p) and pldup = (p,l*inv(d)*u*p).
!
!     for nonsymmetric multicolor Purdue storage format.
!
!  Parameters:
!
!         ndim    row dimension of c,jc arrays
!         jc      integer array giving the column indices of the
!                  corresponding elements in c
!         d       vector of length n giving the diagonal elements
!                  of the matrix
!         c       array of active size n by maxc giving the
!                  off diagonal elements of the matrix.
!                  thus, a = d + c
!         n       order of system
!         ncolor  number of colors used
!         nc      integer vector of length ncolor giving the number
!                  of nodes for each color
!         nt      integer vector of length ncolor giving the number
!                  of upper columns for each color
!         nb      integer vector of length ncolor giving the number
!                  of lower columns for each color
!         p       vector from acceleration algorithm
!         r       workspace vector from acceleration algorithm
!         wksp    workspace vector of length
!                  n + max(nc(i))     if keygs = 1
!                  n                  if keygs = 2
!         pdp     (p,d*p)
!         pldup   (p,l*d*u*p)
!
!  
!
  integer   jc(ndimm,1), nc(1), nt(1), nb(1)
  dimension d(1), c(ndimm,1), p(1), r(1), wksp(1)
!
!  compute pdp = (p,d*p).
!
  n = nn
  ndim = ndimm
  ncolor = ncol
  sum = 0.0
  do 10 i = 1,n
 10   sum = sum + p(i)*d(i)*p(i)
  pdp = sum
!
!  compute pldup = (p,l*inv(d)*u*p) = (u*p,inv(d)*u*p)
!
  np1 = n + 1
  do 15 i = 1,n
 15   r(i) = 0.0
  ist = 1
  do 20 icol = 1,ncolor
     npt = nc(icol)
     mj = nt(icol)
     call vaddp (ndim,ndim,npt,mj,c(ist,1),jc(ist,1),r(ist),p,wksp)
     ist = ist + npt
 20   continue
  do 25 i = 1,n
 25   r(i) = r(i)/d(i)
  do 30 i = 1,n
 30   wksp(i) = 0.0
  ist = 1
  do 35 icol = 1,ncolor
     npt = nc(icol)
     j1 = nt(icol) + 1
     mj = nb(icol)
     call vaddp (ndim,ndim,npt,mj,c(ist,j1),jc(ist,j1),wksp(ist),r,wksp(np1))
     ist = ist + npt
 35   continue
  sum = 0.0
  do 40 i = 1,n
 40   sum = sum + p(i)*wksp(i)
  pldup = sum
  return
end
subroutine sub100 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB100 calls the SRSCP3 preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  call srscp3 (ndim,n,jcoef(ipt1),coef,coef(ipt1),ncolor,iwksp(nc),iwksp(ndt), &
    iwksp(ndb),omega,wksp(iwkpt1),r,z)
  return
end
subroutine sub101 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB101 calls the SRSCP2 preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  call srscp2 (ndim,n,jcoef(ipt1),coef,coef(ipt1),ncolor,iwksp(nc),iwksp(ndt), &
    omega,wksp(iwkpt1),r,z)
  return
end
subroutine sub102 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB102 calls the SRSCP4 preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  call srscp4 (ndim,n,jcoef(ipt1),coef,coef(ipt1),ncolor,iwksp(nc),iwksp(ndt), &
    omega,wksp(iwkpt1),r,z)
  return
end
subroutine sub103 (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SUB103 calls the SSRCP or SSRCPN adaption routine.
!
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
!
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  if (isymm == 0) call ssrcp (ndim,jcoef(ipt1),coef,coef(ipt1), &
     n,ncolor,iwksp(nc),iwksp(ndt),p,r,wksp(iwkpt1),pdp,pldup)
  if (isymm == 1) call ssrcpn (ndim,jcoef(ipt1),coef,coef(ipt1), &
    n,ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),p,r,wksp(iwkpt1),pdp,pldup)
  return
end
subroutine sub104 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB104 calls the ICSCP preconditioner.
!
!     (multicolor Purdue)
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  ipt2 = ifactr + n
  if (propa) call icscp (ndim,ndim,n,jcoef(ipt1),wksp(ifactr), &
    coef(ipt1),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),1,wksp(iwkpt1),r,z)

  if (.not. propa) call icscp (n,ndim,n,jcoef(ipt1),wksp(ifactr), &
    wksp(ipt2),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),0,wksp(iwkpt1),r,z)

  return
end
subroutine sub105 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB105 calls the ICSCPT preconditioner.
!
!     (multicolor Purdue)
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  ipt2 = ifactr + n
  if (propa) call icscpt (ndim,ndim,n,jcoef(ipt1),wksp(ifactr), &
    coef(ipt1),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),1,wksp(iwkpt1),r,z)

  if (.not. propa) call icscpt (n,ndim,n,jcoef(ipt1),wksp(ifactr), &
    wksp(ipt2),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),0,wksp(iwkpt1),r,z)
  return
end
subroutine sub106 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB106 calls the ICSCP1 preconditioner.
!
!     (multicolor Purdue)
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  ipt2 = ifactr + n
  if (propa) call icscp1 (ndim,ndim,n,jcoef(ipt1),wksp(ifactr), &
    coef(ipt1),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),1,wksp(iwkpt1),r,z)

  if (.not. propa) call icscp1 (n,ndim,n,jcoef(ipt1),wksp(ifactr), &
    wksp(ipt2),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),0,wksp(iwkpt1),r,z)

  return
end
subroutine sub107 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB107 calls the ICSCP3 preconditioner.
!
!     (multicolor Purdue)
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  ipt2 = ifactr + n
  if (propa) call icscp3 (ndim,ndim,n,jcoef(ipt1),wksp(ifactr), &
    coef(ipt1),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),1,wksp(iwkpt1),r,z)

  if (.not. propa) call icscp3 (n,ndim,n,jcoef(ipt1),wksp(ifactr), &
    wksp(ipt2),ncolor,iwksp(nc),iwksp(ndt),iwksp(ndb),0,wksp(iwkpt1),r,z)

  return
end
subroutine sub108 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB104 calls the ICSCP2 preconditioner.
!
!     (multicolor Purdue)
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  ipt2 = ifactr + n
  if (propa) call icscp2 (ndim,ndim,n,jcoef(ipt1),wksp(ifactr), &
    coef(ipt1),ncolor,iwksp(nc),iwksp(ndt),1,wksp(iwkpt1),r,z)

  if (.not. propa) call icscp2 (n,ndim,n,jcoef(ipt1),wksp(ifactr), &
    wksp(ipt2),ncolor,iwksp(nc),iwksp(ndt),0,wksp(iwkpt1),r,z)

  return
end
subroutine sub109 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB109 calls the ICSCP4 preconditioner.
!
!     (multicolor Purdue)
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  ipt2 = ifactr + n
  if (propa) call icscp4 (ndim,ndim,n,jcoef(ipt1),wksp(ifactr), &
    coef(ipt1),ncolor,iwksp(nc),iwksp(ndt),1,wksp(iwkpt1),r,z)

  if (.not. propa) call icscp4 (n,ndim,n,jcoef(ipt1),wksp(ifactr), &
    wksp(ipt2),ncolor,iwksp(nc),iwksp(ndt),0,wksp(iwkpt1),r,z)

  return
end
subroutine sub110 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB110 calls PPII, for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba12
!
  call ppii (suba12,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine sub111 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB111 calls PNEU, for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba12
!
  call pneu (suba12,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine sub112 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB112 calls PPII for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba13
!
  call ppii (suba13,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine sub113 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB113 calls PPII for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba14
!
  call ppii (suba14,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine sub114 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB114 calls PNEU for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba13
!
  call pneu (suba13,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine sub115 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUB115 calls PNEU for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba14
!
  call pneu (suba14,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine suba1 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA1 calls MULT2S.
!
  common / dscons / ndim, mdim, maxnz
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mult2s (ndim,maxnz,coef,jcoef,n,x,y)
  return
end
subroutine suba10 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA10 calls RSAP.
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  nr = iwksp(nc)
  nb = iwksp(nc+1)
  nbig = nr + nb
  call rsap (ndim,nbig,n,maxnz,jcoef,coef,x,y,wksp(iwkpt1))
  return
end
subroutine suba11 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA11 calls RSATP.
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  nr = iwksp(nc)
  nb = iwksp(nc+1)
  nbig = nr + nb
  if (isymm == 0) call rsap (ndim,nbig,n,maxnz,jcoef,coef,x,y,wksp(iwkpt1))
  if (isymm == 1) call rsatp (ndim,nbig,n,maxnz,jcoef,coef,x,y,wksp(iwkpt1))
  return
end
subroutine suba12 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA12 calls MULT3.
!
  common / dscons / ndim, mdim, maxnz
  common / cmpart / mpstrt, mpart
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mult3 (mpart,iwksp(mpstrt),coef,jcoef,jcoef(ndim+1),wksp(iwkpt1),x,y)
  return
end
subroutine suba13 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA13 calls MULT3N.
!
  common / dscons / ndim, mdim, maxnz
  common / cmpart / mpstrt, mpart
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mult3n (mpart,iwksp(mpstrt),coef,jcoef,jcoef(ndim+1),wksp(iwkpt1),x,y)
  return
end
subroutine suba14 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA14 calls MUL3NT.
!
  common / dscons / ndim, mdim, maxnz
  common / cmpart / mpstrt, mpart
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mul3nt (mpart,iwksp(mpstrt),coef,jcoef,jcoef(ndim+1),wksp(iwkpt1),x,y)
  return
end
subroutine suba2 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA2 calls MULDC.
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call muldc (ndim,n,coef,ncolor,iwksp(nc),iwksp(maxnew),iwksp(jcnew),x,y)
  return
end
subroutine suba3 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA3 calls MULDCT.
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call muldct (ndim,n,coef,ncolor,iwksp(nc),iwksp(maxnew),iwksp(jcnew),x,y)
  return
end
subroutine suba4 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA4 calls MULT2N.
!
  common / dscons / ndim, mdim, maxnz
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mult2n (ndim,maxnz,coef,jcoef,n,x,y)
  return
end
subroutine suba5 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA5 calls MUL2NT.
!
  common / dscons / ndim, mdim, maxnz
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mul2nt (ndim,maxnz,coef,jcoef,n,x,y)
  return
end
subroutine suba6 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA6 calls RSAD.
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  nr = iwksp(nc)
  nb = iwksp(nc+1)
  nbig = nr + nb
  call rsad (nbig,n,n,ndim,iwksp(maxnew),ndt,ndb,iwksp(jcnew),coef,y,x, &
    wksp(ifactr),wksp(iwkpt1))
  return
end
subroutine suba7 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA7 calls RSATD.
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  nr = iwksp(nc)
  nb = iwksp(nc+1)
  nbig = nr + nb
  call rsatd (nbig,n,n,ndim,iwksp(maxnew),ndt,ndb,iwksp(jcnew),coef,y,x, &
    wksp(ifactr),wksp(iwkpt1))
  return
end
subroutine suba8 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA8 calls MULT1.
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mult1 (ndim,maxnz,coef,jcoef,wksp(iwkpt1),n,x,y)
  return
end
subroutine suba9 (coef,jcoef,wksp,iwksp,n,x,y)
!
!*******************************************************************************
!
!! SUBA9 calls MUL1T.
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension x(1), y(1), coef(1), wksp(1)
!
  call mul1t (ndim,maxnz,coef,jcoef,wksp(iwkpt1),n,x,y)
  return
end
subroutine subq1 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ1 calls PJAC for Jacobi preconditioning.
!
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  call pjac (coef,n,r,z)
  return
end
subroutine subq10 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ10 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxt = maxnz - 1
  call srs2 (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega,irwise, &
    iwksp(iwkpt1),r,z)
  return
end
subroutine subq11 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ11 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxt = maxnz - 1
  call srs4 (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega,irwise, &
    iwksp(iwkpt1),r,z)
  return
end
subroutine subq12 (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SUBQ12 calls the SSOR adaption routine.
!
!
  common / dscons / ndim, mdim, maxnz
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  maxt = maxnz - 1
  call ssord (ndim,maxt,jcoef(2),coef,coef(ndim+1),n,p,r,pdp,pldup)
  return
end
subroutine subq13 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ13 calls ICS for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (propa) call ics (ndim,n,maxt,jcoef(2),wksp(ifactr), &
    coef(ndim+1),1,irwise,iwksp(iwkpt1),r,z)
  if (.not. propa) call ics (n,n,maxt,iwksp(ifacti+1),wksp(ifactr), &
    wksp(ifactr+n),0,irwise,iwksp(iwkpt1),r,z)
  return
end
subroutine subq14 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ14 calls ICS1 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (propa) call ics1 (ndim,n,maxt,jcoef(2),wksp(ifactr), &
    coef(ndim+1),1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call ics1 (n,n,maxt,iwksp(ifacti+1), &
    wksp(ifactr),wksp(ifactr+n), 0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq15 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ15 calls ICS3 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (propa) call ics3 (ndim,n,maxt,jcoef(2),wksp(ifactr), &
    coef(ndim+1),1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call ics3 (n,n,maxt,iwksp(ifacti+1), &
    wksp(ifactr),wksp(ifactr+n),0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq16 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ16 calls ICS2 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (propa) call ics2 (ndim,n,maxt,jcoef(2),wksp(ifactr), &
    coef(ndim+1),1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call ics2 (n,n,maxt,iwksp(ifacti+1), &
    wksp(ifactr),wksp(ifactr+n),0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq17 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ17 calls ICS4 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (propa) call ics4 (ndim,n,maxt,jcoef(2),wksp(ifactr), &
    coef(ndim+1),1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call ics4 (n,n,maxt,iwksp(ifacti+1), &
    wksp(ifactr),wksp(ifactr+n),0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq18 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ18 calls PPII for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba1
!
  call ppii (suba1,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq19 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ19 calls PNEU for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba1
!
  call pneu (suba1,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq2 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ2 calls BDSOL for line Jacobi preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (nstore == 2) isym = 0
  if (nstore == 3) isym = 1
  call bdsol (n,n,kblsz,ndt,ndb,wksp(ifactr),r,z,isym)
  return
end
subroutine subq20 (coef,jcoef,wksp,iwksp,n,u,rhs,unew)
!
!*******************************************************************************
!
!! SUBQ20 calls the basic LSOR iterative step.
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension u(1), rhs(1), unew(1), coef(1), wksp(1)
!
  call sordb (n,ndim,kblsz,kblsz,iwksp(ifacti),lbhb, &
    wksp(ifactr),coef,jcoef,n,omega,u,rhs,unew)

  return
end
subroutine subq21 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ21 calls the LSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim*iwksp(ifacti+2) + 1
  ipt2 = iwksp(ifacti+2) + 1
  call sbsl (n,ndim,n,kblsz,kblsz,lbhb,iwksp(ifacti), &
    wksp(ifactr),coef(ipt1),jcoef(ipt2),r,z,omega,wksp(iwkpt1))

  return
end
subroutine subq22 (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SUBQ22 calls the LSSOR adaption routine.
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  call ssrcd (n,ndim,maxnz,kblsz,iwksp(ifacti),wksp(ifactr), &
    coef,jcoef,n,p,r,wksp(iwkpt1),pdp,pldup)

  return
end
subroutine subq23 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ23 calls PBPII for line LSPOLY preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba1, subq2
!
  call pbpii (suba1,subq2,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg, &
    wksp(iwkpt1),n,r,z)
  return
end
subroutine subq24 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ24 calls PBNEU for line Neumann polynomial preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba1, subq2
!
  call pbneu (suba1,subq2,coef,jcoef,wksp,iwksp,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq25 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ25 calls IBSL for BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  ipt2 = ifactr + n*iwksp(ifacti+2)
  if (lvfill > 0) go to 10
  nwdiag = iwksp(ifacti+2) - ltrunc
  if (propa) call ibsl(n,ndim,n,kblsz,kblsz,lbhb,iwksp(ifacti), &
    wksp(ifactr),coef(ndim*nwdiag+1),jcoef(nwdiag+1),r,z,ivers,wksp(iwkpt1))

  if (.not. propa) call ibsl(n,n,n,kblsz,kblsz,lbhb,iwksp(ifacti), &
    wksp(ifactr),wksp(ipt2),jcoef(nwdiag+1),r,z,ivers,wksp(iwkpt1))

  return
 10   ipt1 = ifacti + 3*lbhb + iwksp(ifacti+2)
  call ibsl (n,n,n,kblsz,kblsz,lbhb,iwksp(ifacti), &
    wksp(ifactr),wksp(ipt2),iwksp(ipt1),r,z,ivers,wksp(iwkpt1))
  return
end
subroutine subq26 (coef,jcoef,wksp,iwksp,n,u,rhs,unew)
!
!*******************************************************************************
!
!! SUBQ26 calls the basic multi-color SOR iterative step
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  integer jcoef(2), iwksp(1)
  dimension u(1), rhs(1), unew(1), coef(1), wksp(1)
!
  call sordmb (n,ndim,n,iwksp(iblock),iwksp(lbhb),ncolor, &
    iwksp(nc),iwksp(ipt),wksp(ifactr),coef,iwksp(jcnew),n,omega,u,rhs,unew)
  return
end
subroutine subq27 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ27 calls the MSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)
  ipt1 = ndim*nwdiag + 1
  ipt2 = ncolor*nwdiag + jcnew
  call sbsln (n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ipt1),iwksp(ipt2),r,z,omega,0,wksp(iwkpt1))
  return
end
subroutine subq28 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ28 calls the MSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)
  ipt1 = ndim*nwdiag + 1
  ipt2 = ncolor*nwdiag + jcnew
  call sbslnt (n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ipt1),iwksp(ipt2),r,z,omega,0,wksp(iwkpt1))
  return
end
subroutine subq29 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ29 calls the MSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)
  ipt1 = ndim*nwdiag + 1
  ipt2 = ncolor*nwdiag + jcnew
  call sbsln1 (n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ipt1),iwksp(ipt2),r,z,omega,0)
  return
end
subroutine subq3 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ3 calls BDSOLT for line Jacobi preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  call bdsolt (n,n,kblsz,ndt,ndb,wksp(ifactr),r,z)
  return
end
subroutine subq30 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ30 calls the MSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)
  ipt1 = ndim*nwdiag + 1
  ipt2 = ncolor*nwdiag + jcnew
  call sbsln3 (n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ipt1),iwksp(ipt2),r,z,omega,0)
  return
end
subroutine subq31 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ31 calls the MSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)
  ipt1 = ndim*nwdiag + 1
  ipt2 = ncolor*nwdiag + jcnew
  call sbsln2 (n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ipt1),iwksp(ipt2),r,z,omega,0,wksp(iwkpt1))
  return
end
subroutine subq32 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ32 calls the MSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)
  ipt1 = ndim*nwdiag + 1
  ipt2 = ncolor*nwdiag + jcnew
  call sbsln4 (n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ipt1),iwksp(ipt2),r,z,omega,0,wksp(iwkpt1))
  return
end
subroutine subq33 (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SUBQ33 calls the MSSOR adaption routine.
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  call ssrcdm (n,ndim,iwksp(lbhb),n,ncolor,iwksp(nc),iwksp(ipt),iwksp(iblock), &
    wksp(ifactr),coef,iwksp(jcnew),n,p,r,wksp(iwkpt1),pdp,pldup)
  return
end
subroutine subq34 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ34 calls IBSLN for multi-color BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)- 2*ltrunc
  if (propa) call ibsln(n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ndim*nwdiag+1),iwksp(jcnew+nwdiag*ncolor), &
    r,z,ivers,0,wksp(iwkpt1))

  if (.not. propa) call ibsln(n,n,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),wksp(iwkpt2), &
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))
  return
end
subroutine subq35 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ35 calls IBSLNT for multi-color bic preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)- 2*ltrunc
  if (propa) call ibslnt(n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),coef(ndim*nwdiag+1),iwksp(jcnew+nwdiag*ncolor), &
    r,z,ivers,0,wksp(iwkpt1))

  if (.not. propa) call ibslnt(n,n,n,n,ncolor,iwksp(nc),iwksp(ipt),iwksp(lbhb), &
    iwksp(iblock),wksp(ifactr),wksp(iwkpt2),iwksp(jcnew+nwdiag*ncolor), &
    r,z,ivers,0,wksp(iwkpt1))
  return
end
subroutine subq36 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ36 calls IBSLN1 for multi-color BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)- 2*ltrunc
  if (propa) call ibsln1(n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),coef(ndim*nwdiag+1),&
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  if (.not. propa) call ibsln1(n,n,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),wksp(iwkpt2), &
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  return
end
subroutine subq37 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ37 calls IBSLN3 for multi-color BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)- 2*ltrunc
  if (propa) call ibsln3(n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),coef(ndim*nwdiag+1), &
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  if (.not. propa) call ibsln3(n,n,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),wksp(iwkpt2), &
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  return
end
subroutine subq38 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ38 calls IBSLN2 for multi-color BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)- 2*ltrunc
  if (propa) call ibsln2(n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt), &
     iwksp(lbhb),iwksp(iblock),wksp(ifactr),coef(ndim*nwdiag+1), &
     iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  if (.not. propa) call ibsln2(n,n,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),wksp(iwkpt2), &
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  return
end
subroutine subq39 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ39 calls IBSLN4 for multi-color BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nwdiag = iwksp(iblock+2) + iwksp(iblock+3*ncolor+2)- 2*ltrunc

  if (propa) call ibsln4(n,ndim,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),coef(ndim*nwdiag+1), &
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  if (.not. propa) call ibsln4(n,n,n,n,ncolor,iwksp(nc),iwksp(ipt), &
    iwksp(lbhb),iwksp(iblock),wksp(ifactr),wksp(iwkpt2), &
    iwksp(jcnew+nwdiag*ncolor),r,z,ivers,0,wksp(iwkpt1))

  return
end
subroutine subq4 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ4 call BMUL or BMULN, for line Jacobi preconditioning
!     (approximate inverse)
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (nstore == 2) isym = 0
  if (nstore == 3) isym = 1
  ift = ifactr + n
  ifb = ifactr + n*(ndt + 1)
  if (isym == 0) call bmul (n,n,ndt,wksp(ifactr),wksp(ift),r,z)
  if (isym == 1) call bmuln (n,n,ndt,ndb,wksp(ifactr),wksp(ift),wksp(ifb),r,z)
  return
end
subroutine subq40 (coef,jcoef,wksp,iwksp,n,u,rhs,unew)
!
!*******************************************************************************
!
!! SUBQ40 calls the basic SOR iterative step
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  integer jcoef(2), iwksp(1)
  dimension u(1), rhs(1), unew(1), coef(1), wksp(1)
!
  maxtp1 = maxt + 1
  call sordn (ndim,n,maxt,maxb,jcoef(2),jcoef(maxt+2),coef,coef(ndim+1), &
    coef(maxtp1*ndim+1),omega,irwise,u,rhs,unew,iwksp(iwkpt1))

  return
end
subroutine subq41 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ41 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxtp1 = maxt + 1
  call srsn (ndim,n,maxt,maxb,jcoef(2),jcoef(maxt+2),coef,coef(ndim+1), &
    coef(ndim*maxtp1+1),omega,irwise,iwksp(iwkpt2),r,z)

  return
end
subroutine subq42 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ42 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxtp1 = maxt + 1
  call srsnt (ndim,n,maxt,maxb,jcoef(2),jcoef(maxt+2),coef,coef(ndim+1), &
    coef(ndim*maxtp1+1),omega,irwise,iwksp(iwkpt2),r,z)

  return
end
subroutine subq43 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ43 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxtp1 = maxt + 1
  call srsn1 (ndim,n,maxb,jcoef(maxt+2),coef,coef(ndim*maxtp1+1),omega,irwise, &
    iwksp(iwkpt2),r,z)
  return
end
subroutine subq44 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ44 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxtp1 = maxt + 1
  call srsn3 (ndim,n,maxb,jcoef(maxt+2),coef,coef(ndim*maxtp1+1),omega,irwise, &
    iwksp(iwkpt2),r,z)
  return
end
subroutine subq45 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ45 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  call srsn2 (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega,irwise, &
    iwksp(iwkpt2),r,z)
  return
end
subroutine subq46 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ46 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  call srsn4 (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega,irwise, &
    iwksp(iwkpt2),r,z)
  return
end
subroutine subq47 (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SUBQ47 calls the SSOR adaption routine.
!
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  maxtp1 = maxt + 1
  call ssordn (ndim,maxt,maxb,jcoef(2),jcoef(maxt+2),coef,coef(ndim+1), &
    coef(ndim*maxtp1+1),n,p,r,wksp(iwkpt1),pdp,pldup)
  return
end
subroutine subq48 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ48 calls ICSN for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  maxtp1 = maxt + 1
  if (propa) call icsn (ndim,n,maxt,maxb,jcoef(2),jcoef(maxt+2), &
    wksp(ifactr),coef(ndim+1),coef(ndim*maxtp1+1),1,irwise, &
    iwksp(iwkpt1),r,z)

  if (.not. propa) call icsn (n,n,maxt,maxb,iwksp(ifacti+1), &
    iwksp(ifacti+maxt+1),wksp(ifactr),wksp(ifactr+n),wksp(ifactr+n*maxtp1), &
    0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq49 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ49 calls ICSNT for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  maxtp1 = maxt + 1
  if (propa) call icsnt (ndim,n,maxt,maxb,jcoef(2),jcoef(maxt+2), &
    wksp(ifactr),coef(ndim+1),coef(ndim*maxtp1+1),1,irwise, &
    iwksp(iwkpt1),r,z)

  if (.not. propa) call icsnt (n,n,maxt,maxb,iwksp(ifacti+1), &
    iwksp(ifacti+maxt+1),wksp(ifactr),wksp(ifactr+n),wksp(ifactr+n*maxtp1), &
    0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq5 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ5 calls BMUL or BMULNT for line Jacobi preconditioning
!     (approximate inverse)
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (nstore == 2) isym = 0
  if (nstore == 3) isym = 1
  ift = ifactr + n
  ifb = ifactr + n*(ndt + 1)
  if (isym == 0) call bmul (n,n,ndt,wksp(ifactr),wksp(ift),r,z)
  if (isym == 1) call bmulnt (n,n,ndt,ndb,wksp(ifactr),wksp(ift),wksp(ifb),r,z)
  return
end
subroutine subq50 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ50 calls ICSN1 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  maxtp1 = maxt + 1
  if (propa) call icsn1 (ndim,n,maxb,jcoef(maxt+2),wksp(ifactr), &
     coef(ndim*maxtp1+1),1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call icsn1 (n,n,maxb,iwksp(ifacti+maxt+1),wksp(ifactr), &
    wksp(ifactr+n*maxtp1),0,irwise,iwksp(iwkpt1),r,z)
  return
end
subroutine subq51 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ51 calls ICSN3 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  maxtp1 = maxt + 1
  if (propa) call icsn3 (ndim,n,maxb,jcoef(maxt+2),wksp(ifactr), &
    coef(ndim*maxtp1+1),1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call icsn3 (n,n,maxb,iwksp(ifacti+maxt+1),wksp(ifactr), &
    wksp(ifactr+n*maxtp1),0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq52 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ52 calls ICSN2 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (propa) call icsn2 (ndim,n,maxt,jcoef(2),wksp(ifactr),coef(ndim+1), &
    1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call icsn2 (n,n,maxt,iwksp(ifacti+1),wksp(ifactr), &
    wksp(ifactr+n),0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq53 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ53 calls ICSN4 for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  if (propa) call icsn4 (ndim,n,maxt,jcoef(2),wksp(ifactr),coef(ndim+1), &
    1,irwise,iwksp(iwkpt1),r,z)

  if (.not. propa) call icsn4 (n,n,maxt,iwksp(ifacti+1),wksp(ifactr), &
    wksp(ifactr+n),0,irwise,iwksp(iwkpt1),r,z)

  return
end
subroutine subq54 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ54 calls PPII for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba4
!
  call ppii (suba4,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq55 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ55 calls PPII for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba5
!
  call ppii (suba5,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq56 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ56 calls PNEU for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba4
!
  call pneu (suba4,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq57 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ57 calls PNEU for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba5
!
  call pneu (suba5,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq58 (coef,jcoef,wksp,iwksp,n,u,rhs,unew)
!
!*******************************************************************************
!
!! SUBQ58 calls the basic LSOR iterative step
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension u(1), rhs(1), unew(1), coef(1), wksp(1)
!
  call sordnb (n,ndim,kblsz,kblsz,iwksp(ifacti),lbhb, &
    wksp(ifactr),coef,jcoef,n,omega,u,rhs,unew)
  return
end
subroutine subq59 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ59 calls the LSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
!
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
  nwdiag = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt1 = ndim*nwdiag + 1
  ipt2 = nwdiag + 1
  call sbsln (n,ndim,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),coef(ipt1),jcoef(ipt2),r,z,omega,1,wksp(iwkpt1))
  return
end
subroutine subq6 (coef,jcoef,wksp,iwksp,n,u,rhs,unew)
!
!*******************************************************************************
!
!! SUBQ6 calls the basic SOR iterative step
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension u(1), rhs(1), unew(1), coef(1), wksp(1)
!
  maxt = maxnz - 1
  call sords (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega, &
    irwise,u,rhs,unew,iwksp(iwkpt1))
  return
end
subroutine subq60 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ60 calls the LSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
!
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
  nwdiag = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt1 = ndim*nwdiag + 1
  ipt2 = nwdiag + 1
  call sbslnt (n,ndim,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),coef(ipt1),jcoef(ipt2),r,z,omega,1,wksp(iwkpt1))
  return
end
subroutine subq61 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ61 calls the LSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
!
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
  nwdiag = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt1 = ndim*nwdiag + 1
  ipt2 = nwdiag + 1
  call sbsln1 (n,ndim,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),coef(ipt1),jcoef(ipt2),r,z,omega,1)
  return
end
subroutine subq62 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ62 calls the LSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
!
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
  nwdiag = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt1 = ndim*nwdiag + 1
  ipt2 = nwdiag + 1
  call sbsln3 (n,ndim,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),coef(ipt1),jcoef(ipt2),r,z,omega,1)
  return
end
subroutine subq63 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ63 calls the LSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
!
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
  nwdiag = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt1 = ndim*nwdiag + 1
  ipt2 = nwdiag + 1
  call sbsln2 (n,ndim,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),coef(ipt1),jcoef(ipt2),r,z,omega,1,wksp(iwkpt1))
  return
end
subroutine subq64 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ64 calls the LSSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
!
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
  nwdiag = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt1 = ndim*nwdiag + 1
  ipt2 = nwdiag + 1
  call sbsln4 (n,ndim,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),coef(ipt1),jcoef(ipt2),r,z,omega,1,wksp(iwkpt1))
  return
end
subroutine subq65 (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SUBQ65 calls the LSSOR adaption routine.
!
!
  common / dscons / ndim, mdim, maxnz
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  call ssrcdn (n,ndim,lbhb,kblsz,iwksp(ifacti),wksp(ifactr),coef,jcoef, &
    n,p,r,wksp(iwkpt1),pdp,pldup)

  return
end
subroutine subq66 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ66 calls PBPII for line LSPOLY preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba4, subq2
!
  call pbpii (suba4,subq2,coef,jcoef,wksp,iwksp,ainf, 0.0,0.0,ndeg, &
    wksp(iwkpt1),n,r,z)
  return
end
subroutine subq67 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ67 calls PBPII for line LSPOLY preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba5, subq3
!
  call pbpii (suba5,subq3,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg, &
    wksp(iwkpt1),n,r,z)
  return
end
subroutine subq68 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ68 calls PBNEU for line Neumann polynomial preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba4, subq2
!
  call pbneu (suba4,subq2,coef,jcoef,wksp,iwksp,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq69 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ69 calls PBNEU for line Neumann polynomial preconditioning.
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba5, subq3
!
  call pbneu (suba5,subq3,coef,jcoef,wksp,iwksp,ndeg,wksp(iwkpt1),n,r,z)
  return
end
subroutine subq7 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ7 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxt = maxnz - 1
  call srs (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega,irwise, &
    iwksp(iwkpt1),r,z)
  return
end
subroutine subq70 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ70 calls IBSLN for BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
!
  n = nn
  nwnew = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt2 = ifactr + n*nwnew
  if (lvfill > 0) go to 10
  nwdiag = nwnew - 2*ltrunc
  if (propa) call ibsln(n,ndim,n,kblsz,1,idumb(1),idumb(2),idumb(3), &
    iwksp(ifacti),wksp(ifactr),coef(ndim*nwdiag+1),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  if (.not. propa) call ibsln(n,n,n,kblsz,1,idumb(1),idumb(2),idumb(3), &
    iwksp(ifacti),wksp(ifactr),wksp(ipt2),jcoef(nwdiag+1),r,z,ivers,1, &
    wksp(iwkpt1))

  return
 10   ipt1 = ifacti + 3*lbhb + nwnew
  call ibsln  (n,n,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),wksp(ipt2),iwksp(ipt1),r,z,ivers,1,wksp(iwkpt1))

  return
end
subroutine subq71 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ71 calls IBSLNT for BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
!
  n = nn
  nwnew = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt2 = ifactr + n*nwnew
  if (lvfill > 0) go to 10
  nwdiag = nwnew - 2*ltrunc
  if (propa) call ibslnt(n,ndim,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),coef(ndim*nwdiag+1),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  if (.not. propa) call ibslnt(n,n,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),wksp(iwkpt2),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  return
 10   ipt1 = ifacti + 3*lbhb + nwnew
  call ibslnt (n,n,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),wksp(ipt2),iwksp(ipt1),r,z,ivers,1,wksp(iwkpt1))

  return
end
subroutine subq72 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ72 calls IBSLN1 for BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
!
  n = nn
  nwnew = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt2 = ifactr + n*nwnew
  if (lvfill > 0) go to 10
  nwdiag = nwnew - 2*ltrunc
  if (propa) call ibsln1(n,ndim,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),coef(ndim*nwdiag+1),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  if (.not. propa) call ibsln1(n,n,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),wksp(iwkpt2),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  return
 10   ipt1 = ifacti + 3*lbhb + nwnew
  call ibsln1 (n,n,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),wksp(ipt2),iwksp(ipt1),r,z,ivers,1,wksp(iwkpt1))

  return
end
subroutine subq73 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ73 calls IBSLN3 for BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
!
  n = nn
  nwnew = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt2 = ifactr + n*nwnew
  if (lvfill > 0) go to 10
  nwdiag = nwnew - 2*ltrunc
  if (propa) call ibsln3(n,ndim,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),coef(ndim*nwdiag+1),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  if (.not. propa) call ibsln3(n,n,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),wksp(iwkpt2),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  return
 10   ipt1 = ifacti + 3*lbhb + nwnew
  call ibsln3 (n,n,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),wksp(ipt2),iwksp(ipt1),r,z,ivers,1,wksp(iwkpt1))

  return
end
subroutine subq74 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ74 calls IBSLN2 for BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
!
  n = nn
  nwnew = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt2 = ifactr + n*nwnew
  if (lvfill > 0) go to 10
  nwdiag = nwnew - 2*ltrunc
  if (propa) call ibsln2(n,ndim,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),coef(ndim*nwdiag+1),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))

  if ( .not. propa ) call ibsln2 (n,n,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),wksp(iwkpt2),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))
  return
 10   ipt1 = ifacti + 3*lbhb + nwnew
  call ibsln2 (n,n,n,kblsz,1,idumb(1),idumb(2),idumb(3),iwksp(ifacti), &
    wksp(ifactr),wksp(ipt2),iwksp(ipt1),r,z,ivers,1,wksp(iwkpt1))
  return
end
subroutine subq75 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ75 calls IBSLN4 for BIC preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  dimension r(1), z(1), coef(1), wksp(1)
  integer idumb(3), jcoef(2), iwksp(1)
  idumb(1) = kblsz
  idumb(2) = 1
  idumb(3) = lbhb
!
  n = nn
  nwnew = iwksp(ifacti+2) + iwksp(ifacti+5)
  ipt2 = ifactr + n*nwnew
  if (lvfill > 0) go to 10
  nwdiag = nwnew - 2*ltrunc
  if (propa) call ibsln4(n,ndim,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),coef(ndim*nwdiag+1),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))
  if (.not. propa) call ibsln4(n,n,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),wksp(iwkpt2),jcoef(nwdiag+1), &
    r,z,ivers,1,wksp(iwkpt1))
  return
 10   ipt1 = ifacti + 3*lbhb + nwnew
  call ibsln4 (n,n,n,kblsz,1,idumb(1),idumb(2), &
    idumb(3),iwksp(ifacti),wksp(ifactr),wksp(ipt2),iwksp(ipt1), &
    r,z,ivers,1,wksp(iwkpt1))
  return
end
subroutine subq76 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ76 calls BDSOL for RS preconditioning.
!
!
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nr = iwksp(nc)
  nb = iwksp(nc+1)
  nbig = nr + nb
  call bdsol (nbig,n,n,ndt,ndb,wksp(ifactr),r,z,1)
  return
end
subroutine subq77 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ77 calls BDSOLT for RS preconditioning.
!
!
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  n = nn
  nr = iwksp(nc)
  nb = iwksp(nc+1)
  nbig = nr + nb
  call bdsolt (nbig,n,n,ndt,ndb,wksp(ifactr),r,z)
  return
end
subroutine subq78 (coef,jcoef,wksp,iwksp,n,u,rhs,unew)
!
!*******************************************************************************
!
!! SUBQ78 calls the basic SOR iterative step
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  integer jcoef(2), iwksp(1)
  dimension u(1), rhs(1), unew(1), coef(1), wksp(1)
!
  ip1 = ndim + 1
  ip2 = ndim*(maxt + 1) + 1
  call sorp (ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2),coef, &
    coef(ip1),coef(ip2),omega,u,rhs,unew)
  return
end
subroutine subq79 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ79 calls the SRSP preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ip1 = ndim + 1
  ip2 = ndim*(maxt + 1) + 1
  call srsp (ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2),coef,coef(ip1), &
    coef(ip2),omega,r,z)
  return
end
subroutine subq8 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ8 calls the SRS1 preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxt = maxnz - 1
  call srs1 (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega,irwise, &
    iwksp(iwkpt1),r,z)
  return
end
subroutine subq80 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ80 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ip1 = ndim + 1
  ip2 = ndim*(maxt + 1) + 1
  call srsntp (ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2),coef,coef(ip1), &
    coef(ip2),omega,r,z)
  return
end
subroutine subq81 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ81 calls the SRSP1 preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ip2 = ndim*(maxt + 1) + 1
  call srsp1 (ndim,n,maxb,jcoef(ip2),coef,coef(ip2),omega,r,z)
  return
end
subroutine subq82 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ82 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ip2 = ndim*(maxt + 1) + 1
  call srsp3 (ndim,n,maxb,jcoef(ip2),coef,coef(ip2),omega,r,z)
  return
end
subroutine subq83 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ83 calls the SRSP2 preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ip1 = ndim + 1
  call srsp2 (ndim,n,maxt,jcoef(ip1),coef,coef(ip1),omega,r,z)
  return
end
subroutine subq84 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ84 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ip1 = ndim + 1
  call srsp4 (ndim,n,maxt,jcoef(ip1),coef,coef(ip1),omega,r,z)
  return
end
subroutine subq85 (coef,jcoef,wksp,iwksp,n,p,r,pdp,pldup)
!
!*******************************************************************************
!
!! SUBQ85 calls the SSOR adaption routine.
!
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension p(1), r(1), coef(1), wksp(1)
!
  ip1 = ndim + 1
  ip2 = ndim*(maxt + 1) + 1
  if (isymm == 0) call ssorp (ndim,maxt,jcoef(ip1),coef, &
    coef(ip1),n,p,r,wksp(iwkpt1),pdp,pldup)
  if (isymm /= 0) call ssorpn (ndim,maxt,maxb,jcoef(ip1), &
    jcoef(ip2),coef,coef(ip1),coef(ip2),n,p,r,wksp(iwkpt1),pdp,pldup)
  return
end
subroutine subq86 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ86 calls ICS for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  logical symm
!
  n = nn
  symm = isymm == 0
  if (.not. propa) go to 10
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     if (symm) call icsp (ndim,ndim,n,maxt,jcoef(ip1), &
       wksp(ifactr),coef(ip1),1,r,z)

     if (.not. symm) call icsnp (ndim,ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2), &
       wksp(ifactr),coef(ip1),coef(ip2),1,r,z)
     return
 10   if (lvfill > 0) go to 15
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsnp (n,ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2), &
       wksp(ifactr),wksp(ip3),wksp(ip4),0,r,z)
     return
 15   continue
     ip1 = ifacti + n
     ip2 = ifacti + n*(maxt + 1)
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp (n,n,n,maxt,iwksp(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsnp (n,n,n,maxt,maxb,iwksp(ip1),iwksp(ip2), &
       wksp(ifactr),wksp(ip3),wksp(ip4),0,r,z)
     return
end
subroutine subq87 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ87 calls ICS for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  logical symm
!
  n = nn
  symm = isymm == 0
  if (.not. propa) go to 10
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     if (symm) call icsp (ndim,ndim,n,maxt,jcoef(ip1),wksp(ifactr),coef(ip1), &
       1,r,z)

     if (.not. symm) call icsntp (ndim,ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2), &
       wksp(ifactr),coef(ip1),coef(ip2),1,r,z)

     return
 10   if (lvfill > 0) go to 15
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr),wksp(ip3),0,r,z)

     if (.not. symm) call icsntp (n,ndim,n,maxt,maxb,jcoef(ip1),jcoef(ip2), &
       wksp(ifactr),wksp(ip3),wksp(ip4),0,r,z)

     return
 15   continue
     ip1 = ifacti + n
     ip2 = ifacti + n*(maxt + 1)
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp (n,n,n,maxt,iwksp(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsntp (n,n,n,maxt,maxb,iwksp(ip1),iwksp(ip2), &
       wksp(ifactr),wksp(ip3),wksp(ip4),0,r,z)
     return
end
subroutine subq88 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ88 calls ICS for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  logical symm
!
  n = nn
  symm = isymm == 0
  if (.not. propa) go to 10
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     if (symm) call icsp1 (ndim,ndim,n,maxt,jcoef(ip1),wksp(ifactr), &
       coef(ip1),1,r,z)

     if (.not. symm) call icsnp1 (ndim,ndim,n,maxb,jcoef(ip2),wksp(ifactr), &
       coef(ip2),1,r,z)

     return
 10   if (lvfill > 0) go to 15
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp1 (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr),wksp(ip3), &
       0,r,z)
     if (.not. symm) call icsnp1 (n,ndim,n,maxb,jcoef(ip2),wksp(ifactr), &
       wksp(ip4),0,r,z)
     return
 15   continue
     ip1 = ifacti + n
     ip2 = ifacti + n*(maxt + 1)
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp1 (n,n,n,maxt,iwksp(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsnp1 (n,n,n,maxb,iwksp(ip2),wksp(ifactr), &
       wksp(ip4),0,r,z)
  return
end
subroutine subq89 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ89 calls ICS for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  logical symm
!
  n = nn
  symm = isymm == 0
  if (.not. propa) go to 10
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     if (symm) call icsp3 (ndim,ndim,n,maxt,jcoef(ip1),wksp(ifactr),coef(ip1), &
       1,r,z)
     if (.not. symm) call icsnp3 (ndim,ndim,n,maxb,jcoef(ip2),wksp(ifactr), &
       coef(ip2),1,r,z)
     return
 10   if (lvfill > 0) go to 15
     ip1 = ndim + 1
     ip2 = ndim*(maxt + 1) + 1
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp3 (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr),wksp(ip3), &
       0,r,z)
     if (.not. symm) call icsnp3 (n,ndim,n,maxb,jcoef(ip2),wksp(ifactr), &
       wksp(ip4),0,r,z)
     return
 15   continue
     ip1 = ifacti + n
     ip2 = ifacti + n*(maxt + 1)
     ip3 = ifactr + n
     ip4 = n*(maxt + 1)+ ifactr
     if (symm) call icsp3 (n,n,n,maxt,iwksp(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsnp3 (n,n,n,maxb,iwksp(ip2),wksp(ifactr), &
       wksp(ip4),0,r,z)
     return
end
subroutine subq9 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ9 calls the SSOR preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  maxt = maxnz - 1
  call srs3 (ndim,n,maxt,jcoef(2),coef,coef(ndim+1),omega,irwise, &
    iwksp(iwkpt1),r,z)
  return
end
subroutine subq90 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ90 calls ICS for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  logical symm
!
  n = nn
  symm = isymm == 0
  if (.not. propa) go to 10
     ip1 = ndim + 1
     if (symm) call icsp2 (ndim,ndim,n,maxt,jcoef(ip1),wksp(ifactr),coef(ip1), &
       1,r,z)
     if (.not. symm) call icsnp2 (ndim,ndim,n,maxt,jcoef(ip1),wksp(ifactr), &
       coef(ip1),1,r,z)
     return
 10   if (lvfill > 0) go to 15
     ip1 = ndim + 1
     ip3 = ifactr + n
     if (symm) call icsp2 (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsnp2 (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr), &
       wksp(ip3),0,r,z)
     return
 15   continue
     ip1 = ifacti + n
     ip3 = ifactr + n
     if (symm) call icsp2 (n,n,n,maxt,iwksp(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsnp2 (n,n,n,maxt,iwksp(ip1),wksp(ifactr), &
       wksp(ip3),0,r,z)
     return
end
subroutine subq91 (coef,jcoef,wksp,iwksp,nn,r,z)
!
!*******************************************************************************
!
!! SUBQ91 calls ICS for IC(S) preconditioning.
!
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  common / cfactr / nfactr, nfacti, ifactr, ifacti, timfac
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  logical symm
!
  n = nn
  symm = isymm == 0
  if (.not. propa) go to 10
     ip1 = ndim + 1
     if (symm) call icsp4 (ndim,ndim,n,maxt,jcoef(ip1), &
       wksp(ifactr),coef(ip1),1,r,z)

     if (.not. symm) call icsnp4 (ndim,ndim,n,maxt,jcoef(ip1),wksp(ifactr), &
       coef(ip1),1,r,z)
     return
 10   if (lvfill > 0) go to 15
     ip1 = ndim + 1
     ip3 = ifactr + n
     if (symm) call icsp4 (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr),wksp(ip3), &
       0,r,z)
     if (.not. symm) call icsnp4 (n,ndim,n,maxt,jcoef(ip1),wksp(ifactr), &
       wksp(ip3),0,r,z)
     return
 15   continue
     ip1 = ifacti + n
     ip3 = ifactr + n
     if (symm) call icsp4 (n,n,n,maxt,iwksp(ip1),wksp(ifactr),wksp(ip3),0,r,z)
     if (.not. symm) call icsnp4 (n,n,n,maxt,iwksp(ip1),wksp(ifactr), &
       wksp(ip3),0,r,z)
     return
end
subroutine subq92 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ92 calls PPII for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba8
!
  call ppii (suba8,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine subq93 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ93 calls PPII for LSPOLY preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  common / itcom8 / ainf
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba9
!
  call ppii (suba9,coef,jcoef,wksp,iwksp,ainf,0.0,0.0,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine subq94 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ94 calls PNEU for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba8
!
  call pneu (suba8,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine subq95 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ95 calls PNEU for Neumann polynomial preconditioning.
!
!
  common / itcom6 / method, iscale, iperm, nstore, ifact, kblsz, lvfill, &
         ltrunc, ndeg, ipropa, isymm, ifctv
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
  external suba9
!
  call pneu (suba9,coef,jcoef,wksp,iwksp,coef,ndeg,wksp(iwkpt2),n,r,z)
  return
end
subroutine subq96 (coef,jcoef,wksp,iwksp,n,u,rhs,unew)
!
!*******************************************************************************
!
!! SUBQ96 calls the basic multi-color SOR iterative step
!
  common / dscons / ndim, mdim, maxnz
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
  integer jcoef(2), iwksp(1)
  dimension u(1), rhs(1), unew(1), coef(1), wksp(1)
!
  call sorcp (ndim,n,jcoef(ndim+1),coef,coef(ndim+1),ncolor, &
    iwksp(nc),iwksp(ndt),iwksp(ndb),omega,u,rhs,unew)
  return
end
subroutine subq97 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ97 calls the SRSCP preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  call srscp (ndim,n,jcoef(ipt1),coef,coef(ipt1),ncolor,iwksp(nc), &
    iwksp(ndt),iwksp(ndb),omega,wksp(iwkpt1),r,z)
  return
end
subroutine subq98 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ98 calls the SRSCPT preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  call srscpt (ndim,n,jcoef(ipt1),coef,coef(ipt1),ncolor,iwksp(nc),iwksp(ndt), &
    iwksp(ndb),omega,wksp(iwkpt1),r,z)
  return
end
subroutine subq99 (coef,jcoef,wksp,iwksp,n,r,z)
!
!*******************************************************************************
!
!! SUBQ99 calls the SRSCP1 preconditioner.
!
!
! *** begin -- itpack common
!
  logical           omgadp
  common / itcom5 / omgadp, omega, alphab, betab, fff, specr
!
! *** end   -- itpack common
!
  common / dscons / ndim, mdim, maxnz
  common / point  / iptscl, iwkpt1, iwkpt2, iwkpt3
  logical           propa
  common / cblock / propa, ncolor, maxd, nc, ipt, maxnew, jcnew, lbhb, iblk, &
         ncmax
  common / intern / ndt, ndb, maxt, maxb, ivers, irwise
  integer jcoef(2), iwksp(1)
  dimension r(1), z(1), coef(1), wksp(1)
!
  ipt1 = ndim + 1
  call srscp1 (ndim,n,jcoef(ipt1),coef,coef(ipt1),ncolor,iwksp(nc),iwksp(ndt), &
    iwksp(ndb),omega,wksp(iwkpt1),r,z)
  return
end
subroutine t1prod (lda,ldb,ldc,ldd,ldj,nn,np,nq,ma,mb,mc,md,incb,incc,incd, &
  ja,jb,jc,jd,a,b,c,d)
!
!*******************************************************************************
!
!! T1PROD computes D = D - C*A*B restricted to the sparsity pattern of D.  
!
! a is assumed to
!  be in nonsymmetric storage mode.
!
!      c is np x nn        b is nn x nq
!      a is nn x nn        d is np x nq
!
!  definition of parameters --
!
!         lda,ldb,      row dimension of arrays a,b,c,d
!          ldc,ldd
!         ldj           row dimension of arrays ja,jb,jc,jd
!         nn,np,nq      orders of arrays
!         ma,mb,mc,md   columns (diagonals) in arrays a,b,c,d
!         incb,incc,    offsets for diagonal numbers of b,c,d arrays
!           incd
!         ja,jb,jc,jd   diagonal index arrays for a,b,c,d
!         a,b,c,d       arrays of dimension n x (ma,mb,mc, or md)
!
!  
!
  integer ja(ldj,1), jb(ldj,1), jc(ldj,1), jd(ldj,1)
  dimension a(lda,1), b(ldb,1), c(ldc,1), d(ldd,1)
!
  n = nn
  do 40 lc = 1,mc
     i = jc(1,lc) - incc
     ia1 = max (1,1-i)
     ib1 = min (np,n-i)
     do 35 la = 1,ma
        j = ja(1,la)
        l1 = i + j
        ia2 = max (ia1,1-l1)
        ib2 = min (ib1,n-l1)
        do 30 lb = 1,mb
           k = jb(1,lb) - incb
           l = l1 + k
           do 15 ld = 1,md
              if (jd(1,ld)-incd == l) go to 20
 15            continue
           go to 30
 20            ist = max (ia2,1-l)
           ied = min (ib2,nq-l)
           do 25 m = ist,ied
 25            d(m,ld) = d(m,ld) - c(m,lc)*a(m+i,la)*b(m+l1,lb)
 30         continue
 35      continue
 40   continue
  return
end
subroutine t2prod (nn,nda,ndb,ndc,ndd,ma,mbb,mc,mdd,incb,incc,incd,ja,jb,jc, &
  jd,a,b,c,d)
!
!*******************************************************************************
!
!! T2PROD computes D = D - (C**t)*A*B restricted to the sparsity pattern of D.  
!
!  a is assumed to be symmetric.
!
!  Parameters:
!
!         n             orders of arrays a,b,c,d
!         nda,ndb,ndc,  row dimensions of arrays a,b,c,d
!          ndd
!         ma,mb,mc,md   columns (diagonals) in arrays a,b,c,d
!         incb,incc,    offsets for diagonal numbers of b,c,d arrays
!           incd
!         ja,jb,jc,jd   diagonal index arrays for a,b,c,d
!         a,b,c,d       arrays of dimension n x (ma,mb,mc, or md)
!
!  
!
  integer   ja(1), jb(1), jc(1), jd(1)
  dimension a(nda,1), b(ndb,1), c(ndc,1), d(ndd,1)
!
  n = nn
  mb = mbb
  md = mdd
  do 65 lc = 1,mc
     i = jc(lc) - incc
     ia1 = max (1,i+1)
     ib1 = min (n,n+i)
     do 60 la = 1,ma
        j = ja(la)
        l1 = -i + j
        ia2 = max (ia1,1-l1)
        ib2 = min (ib1,n-l1)
        do 30 lb = 1,mb
           k = jb(lb) - incb
           l = l1 + k
           do 15 ld = 1,md
              if (jd(ld)-incd == l) go to 20
 15            continue
           go to 30
 20            ist = max (ia2,1-l)
           ied = min (ib2,n-l)
           do 25 ir = ist,ied
 25            d(ir,ld) = d(ir,ld) - c(ir-i,lc)*a(ir-i,la)*b(ir+l1,lb)
 30         continue
        if (j == 0) go to 60
        l1 = -i - j
        ia2 = max (ia1,1-l1)
        ib2 = min (ib1,n-l1)
        do 55 lb = 1,mb
           k = jb(lb) - incb
           l = l1 + k
           do 40 ld = 1,md
              if (jd(ld)-incd == l) go to 45
 40            continue
           go to 55
 45            ist = max (ia2,1-l)
           ied = min (ib2,n-l)
           do 50 ir = ist,ied
 50            d(ir,ld) = d(ir,ld) - c(ir-i,lc)*a(ir+l1,la)*b(ir+l1,lb)
 55         continue
 60      continue
 65   continue
  return
end
function tau (ii)
!
!*******************************************************************************
!
!! TAU sets TAU for the SOR method.
!
!
!  Parameters:
!
!          ii     number of times parameters have been changed
!
!  
!
!
  dimension t(9)
!
  data  t(1), t(2), t(3), t(4), t(5), t(6),  t(7),  t(8), t(9) &
    / 1.5, 1.8, 1.85, 1.9, 1.94, 1.96, 1.975, 1.985, 1.992 /
!
  tau = t(9)
  if (ii <= 8) tau = t(ii)
  return
end
subroutine tbs (n,t,x)
!
!*******************************************************************************
!
!! TBS does a back substitution.
!
!
!  This has the form  (i + t)*x = y  where t is the
!     first super-diagonal.
!
!  Parameters:
!
!          n      order of the system
!          t      vector of length n-1 containing the super-
!                  diagonal elements
!          x      on input, x contains y
!                 on output, x contains the solution to (i - t)*x = y
!
!  
!
  dimension t(1), x(1)
!
  do 10 i = n-1,1,-1
 10   x(i) = x(i) - t(i)*x(i+1)
  return
end
subroutine tbsm (nn,nsize,t,x)
!
!*******************************************************************************
!
!! TBSM does a back substitution.
!
!  This has the form  (i + t)*x = y  where t
!     is a super diagonal composed of independent subsystems of
!     size nsize.
!
!  Parameters:
!
!          n      order of system
!          nsize  order of the individual subsystems
!          t      linear array of length n-1 containing the super-
!                  diagonal elements of the factorizations
!          x      on input, x contains y
!                 the solution to (i + t)*x = y
!
!  
!
  dimension t(nsize,1), x(nsize,1)
!
  n = nn
  nsys = n/nsize
  do 15 i = nsize-1,1,-1
     do 10 j = 1,nsys
 10      x(i,j) = x(i,j) - t(i,j)*x(i+1,j)
 15   continue
  return
end
subroutine tfac (nn,d,t)
!
!*******************************************************************************
!
!! TFAC computes a factorization of a single symmetric tridiagonal matrix.
!
!
!  The matrix is contained in d and t and the factorization overwrites it.
!
!  Parameters:
!
!          n      order of system (= nn)
!          d      vector of length n containing the diagonal
!                  elements of the matrix
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the matrix
!
!  
!
  dimension d(1), t(1)
!
  n = nn
  nm1 = n - 1
  do 10 i = 2,n
 10   d(i) = d(i) - (t(i-1)*t(i-1))/d(i-1)
  do 15 i = 1,n
 15   d(i) = 1.0/d(i)
  do 20 i = 1,nm1
 20   t(i) = d(i)*t(i)
  return
end
subroutine tfacm (nn,nsize,d,t)
!
!*******************************************************************************
!
!! TFACM factors of multiple independent symmetric tridiagonal matrices.
!
!
!  The matrices are contained in d and t.
!
!  Parameters:
!
!          n      order of global system (= nn)
!          nsize  size of the individual subsystems
!          d      linear array of length n containing the
!                  diagonal elements of the systems
!          t      linear array of length n-1 containing the
!                  super-diagonal elements of the systems
!
!  
!
  dimension d(nsize,1), t(nsize,1)
!
  n = nn
  nm1 = n - 1
  nsys = n/nsize
  do 10 i = 2,nsize
     do 5 j = 1,nsys
 5       d(i,j) = d(i,j) - (t(i-1,j)**2)/d(i-1,j)
 10   continue
  call vinv (n,d)
  call vexopy (nm1,t,d,t,3)
  return
end
subroutine tfacn (nn,d,t,b)
!
!*******************************************************************************
!
!! TFACN factors a nonsymmetric tridiagonal matrix.
!
!
!  The matrix is contained in d, t, and b and the factorization
!     replaces it.
!
!  Parameters:
!
!          n      order of system (= nn)
!          d      vector of length n containing the diagonal
!                  elements of the matrix
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the matrix
!          b      vector of length n-1 containing the sub-
!                  diagonal elements of the matrix
!
!  
!
  dimension d(1), t(1), b(1)
!
  n = nn
  nm1 = n - 1
  do 10 i = 2,n
 10   d(i) = d(i) - b(i-1)*t(i-1)/d(i-1)
  do 15 i = 1,n
 15   d(i) = 1.0/d(i)
  do 20 i = 1,nm1
     t(i) = d(i)*t(i)
     b(i) = d(i)*b(i)
 20   continue
  return
end
subroutine tfacnm (nn,nsize,d,t,b)
!
!*******************************************************************************
!
!! TFACNM factors multiple independent nonsymmetric tridiagonal matrices.
!
!
!  The matrices are contained in
!     d, t, and b.
!
!  Parameters:
!
!          n      order of global system (= nn)
!          nsize  order of single subsystem
!          d      linear array of length n containing the
!                  diagonal elements of the systems
!          t      linear array of length n-1 containing the
!                  super-diagonal elements of the systems
!          b      linear array of length n-1 containing the
!                  sub-diagonal elements of the systems
!
!  
!
  dimension d(nsize,1), t(nsize,1), b(nsize,1)
!
  n = nn
  nm1 = n - 1
  nsys = n/nsize
  do 10 i = 2,nsize
     do 5 j = 1,nsys
 5       d(i,j) = d(i,j) - b(i-1,j)*t(i-1,j)/d(i-1,j)
 10   continue
  call vinv (n,d)
  call vexopy (nm1,t,d,t,3)
  call vexopy (nm1,b,d,b,3)
  return
end
subroutine tfs (n,b,x)
!
!*******************************************************************************
!
!! TFS does a forward substitution.
!
!
!  This has the form (i + b)*x = y,
!     where b is the first sub-diagonal.
!
!  Parameters:
!
!          n      order of system
!          b      vector of length n-1 containing the sub-
!                  diagonal elements
!          x      on input, x contains y
!                 on output, x contains the solution to (i - b)*x = y
!
!  
!
  dimension b(1), x(1)
!
  do 10 i = 2,n
 10   x(i) = x(i) - b(i-1)*x(i-1)
  return
end
subroutine tfsm (nn,nsize,b,x)
!
!*******************************************************************************
!
!! TFSM does a forward substitution.
!
!
!  This has the form  (i + b)*x = y  where b
!     is a sub-diagonal composed of independent subsystems of
!     size nsize.
!
!  Parameters:
!
!          n      order of system
!          nsize  order of the individual subsystems
!          b      linear array of length n-1 containing the sub-
!                  diagonal elements of the factorizations
!          x      on input, x contains y
!                 on output, x contains the solution to (i + b)*x = y
!
!  
!
  dimension b(nsize,1), x(nsize,1)
!
  n = nn
  nsys = n/nsize
  do 20 i = 2,nsize
     do 15 j = 1,nsys
 15      x(i,j) = x(i,j) - b(i-1,j)*x(i-1,j)
 20   continue
  return
end
function timer(timdmy)
!
!*******************************************************************************
!
!! TIMER is a routine to return the execution time in seconds.  
!
!
  real tarray(2)
  real timdmy
  real timer
!
  call etime(tarray)
  timdmy=tarray(1)+tarray(2)
  timer=timdmy

  return
end
subroutine tinv (nn,d,t)
!
!*******************************************************************************
!
!! TINV computes an approximate inverse to a single tridiagonal symmetric matrix.  
!
!
!  d and u must contain upon input the
!     output from a factorization routine.
!
!  Parameters:
!
!          n      order of system (= nn)
!          d      vector of length n containing the diagonal
!                  elements of the factorization
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the factorization
!
!  
!
  dimension d(1), t(1)
!
  n = nn
  nm1 = n - 1
!
  do 10 i = nm1,1,-1
 10   d(i) = d(i) + t(i)*t(i)*d(i+1)
  do 15 i = 1,nm1
 15   t(i) = -d(i+1)*t(i)
  return
end
subroutine tinvm (nn,nsize,d,t)
!
!*******************************************************************************
!
!! TINVM computes an approximate inverse to multiple tridiagonal symmetric matrices.  
!
!  d and t must contain upon input the
!     output from a factorization routine.
!
!  Parameters:
!
!          n      order of system (= nn)
!          nsize  size of a single subsystem
!          d      vector of length n containing the diagonal
!                  elements of the factorization
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the factorization
!
!  
!
  dimension d(nsize,1), t(nsize,1)
!
  n = nn
  nm1 = n - 1
  nsys = n/nsize
  nsm1 = nsize - 1
!
  do 20 i = nsm1,1,-1
     do 15 l = 1,nsys
 15      d(i,l) = d(i,l) + t(i,l)*t(i,l)*d(i+1,l)
 20   continue
  call vemxty (nm1,t,d(2,1),t)
  return
end
subroutine tinvn (nn,d,t,b)
!
!*******************************************************************************
!
!! TINVN computes an approximate inverse to a single tridiagonal nonsymmetric matrix.  
!
!
!  d, b, and t must contain upon
!     input the output from a factorization routine.
!
!  Parameters:
!
!          n      order of system (= nn)
!          d      vector of length n containing the diagonal
!                  elements of the factorization
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the factorization
!          b      vector of length n-1 containing the sub-
!                  diagonal elements of the factorization
!
!  
!
  dimension d(1), t(1), b(1)
!
  n = nn
  nm1 = n - 1
!
  do 10 i = nm1,1,-1
 10   d(i) = d(i) + b(i)*t(i)*d(i+1)
  do 20 i = 1,nm1
     t(i) = -d(i+1)*t(i)
     b(i) = -d(i+1)*b(i)
 20   continue
  return
end
subroutine tinvnm (nn,nsize,d,t,b)
!
!*******************************************************************************
!
!! TINVNM computes an approximate inverse to multiple tridiagonal nonsymmetric matrices.  
!
!
!  d, t, and b must contain upon
!     input the output from a factorization routine.
!
!  Parameters:
!
!          n      order of system (= nn)
!          nsize  size of a single subsystem
!          d      vector of length n containing the diagonal
!                  elements of the factorization
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the factorization
!          b      vector of length n-1 containing the sub-
!                  diagonal elements of the factorization
!
!  
!
  dimension d(nsize,1), t(nsize,1), b(nsize,1)
!
  n = nn
  nm1 = n - 1
  nsys = n/nsize
  nsm1 = nsize - 1
!
  do 20 i = nsm1,1,-1
     do 15 l = 1,nsys
 15      d(i,l) = d(i,l) + b(i,l)*t(i,l)*d(i+1,l)
 20   continue
  call vemxty (nm1,t,d(2,1),t)
  call vemxty (nm1,b,d(2,1),b)
  return
end
subroutine tmult (n,nblk,nband,ctac,eta,u,v)
!
!*******************************************************************************
!
!! TMULT omputes the product of the t-matrix with a vector.  
!
!
!  here, t = c*((c**t)*a*c)**(-1) * c**t, a projection.
!
  dimension ctac(nblk,nband), eta(1), u(1), v(1)
!
  nbsiz = n / nblk
  nhband = (nband-1)/2
  nhbp1 = nhband + 1
!
! form the eta vector - aggregation step.
!
  do 1 i=0,nblk-1
!1    eta(i) = vadd (nbsiz,u(1+i*nbsiz))
  ip1 = i + 1
  eta(ip1) = 0e0
  do 1 j=1,nbsiz
 1    eta(ip1) = eta(ip1) + u(i*nbsiz+j)
!
! perform the forward solve.
!
  if (nhband == 0) go to 40
  do 2 irow=2,nblk
  ibeg = max (1,irow-nhband)
  iend = irow - 1
  ind = nhbp1 - irow
  do 3 icol = ibeg,iend
 3    eta(irow) = eta(irow) - eta(icol)*ctac(irow,ind+icol)
 2    continue
!
! perform the diagonal solve.
!
 40   do 4 i=1,nblk
 4    eta(i) = eta(i) / ctac(i,nhbp1)
!
! perform the back solve.
!
  if (nhband == 0) go to 41
  do 5 i=1,nblk-1
  irow = nblk - i
  ibeg = irow + 1
  iend = min (irow+nhband,nblk)
  ind = nhbp1 - irow
  do 6 icol = ibeg,iend
 6    eta(irow) = eta(irow) - eta(icol)*ctac(irow,ind+icol)
 5    continue
!
! form the vector t*u - disaggregation step.
!
 41   do 7 i=0,nblk-1
  val = eta(i+1)
!7    call vfill (nbsiz,v(1+i*nbsiz),eta(i+1))
  do 7 j=1,nbsiz
 7    v(i*nbsiz+j) = val
!
  return
end
subroutine tsoln (nn,d,t,b,y,x)
!
!*******************************************************************************
!
!! TSOLN solves A*x = y for x, for a tridiagonal system.  
!
!
!  d, t, and b contain
!     the main diagonal, the first super-diagonal, and the first
!     sub-diagonal, respectively of the factorization.
!
!  Parameters:
!
!          n      order of system
!          d      vector of length n containing the diagonal
!                  elements of the factorization matrix
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the factorization
!          b      vector of length n-1 containing the sub-
!                  diagonal elements of the factorization
!          y      the right-hand side
!          x      the solution to ax = y
!
!  
!
  dimension d(1), t(1), b(1), x(1), y(1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call tfs (n,b,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call tbs (n,t,x)
  return
end
subroutine tsolnm (nn,nsize,d,t,b,y,x)
!
!*******************************************************************************
!
!! TSOLNM solves A*x = y for x, for multiple tridiagonal systems A.  
!
!
!  d, t, and b contain
!     the main diagonal, the first super-diagonal, and the first
!     sub-diagonal, respectively of the factorization.
!
!  Parameters:
!
!          n      order of system
!          nsize  size of an individual subsystem
!          d      vector of length n containing the diagonal
!                  elements of the factorization matrix
!          t      vector of length n-1 containing the super-
!                  diagonal elements of the factorization
!          b      vector of length n-1 containing the sub-
!                  diagonal elements of the factorization
!          y      the right-hand side
!          x      the solution to ax = y
!
!  
!
  dimension d(1), t(1), b(1), x(1), y(1)
!
  n = nn
  do 10 i = 1,n
 10   x(i) = y(i)
  call tfsm (n,nsize,b,x)
  do 15 i = 1,n
 15   x(i) = d(i)*x(i)
  call tbsm (n,nsize,t,x)
  return
end
subroutine tsum (nn,lda,ldb,ldc,ma,mbb,mc,mdd,incb,incc,incdd,ja,jb,jc,jd, &
  a,b,c,rows,cols,wksp,icodee,omegaa)
!
!*******************************************************************************
!
!! TSUM computes the row and column sum of (C**t)*A*B restricted to the sparsity pattern of JD.  
!
!
!  a is assumed to be symmetric.
!
!  Parameters:
!
!         n             orders of arrays a,b,c,d
!         lda,ldb,ldc   row dimensions of arrays a,b,c
!         ma,mb,mc,md   columns (diagonals) in arrays a,b,c,d
!         incb,incc,    offsets for diagonal numbers of b,c,d arrays
!           incd
!         ja,jb,jc,jd   diagonal index arrays for a,b,c,d
!         a,b,c         arrays of dimension n x (ma,mb,md)
!         rows          row sum of d = (c**t)*a*b upon output
!         cols          column sum of d upon output
!         wksp          workspace array of length n
!         icode         key
!                        = 0  if c /= b
!                        = 1  if c == b
!         omega         relaxation factor between 0 and 1
!
!  
!
  integer   ja(1), jb(1), jc(1), jd(1)
  dimension a(lda,1), b(ldb,1), c(ldc,1), wksp(1), rows(1), cols(1)
!
  n = nn
  mb = mbb
  md = mdd
  incd = incdd
  icode = icodee
  omega = omegaa
  do 95 lc = 1,mc
     i = jc(lc) - incc
     ia1 = max (1,i+1)
     ib1 = min (n,n+i)
     do 90 la = 1,ma
        j = ja(la)
        l1 = -i + j
        ia2 = max (ia1,1-l1)
        ib2 = min (ib1,n-l1)
        do 45 lb = 1,mb
           k = jb(lb) - incb
           l = l1 + k
           do 10 ld = 1,md
              if (jd(ld)-incd == l) go to 15
 10            continue
           go to 45
 15            ist = max (ia2,1-l)
           ied = min (ib2,n-l)
           do 20 kk = ist,ied
 20            wksp(kk-ist+1) = c(kk-i,lc)*a(kk-i,la)*b(kk+l1,lb)
           do 25 kk = ist,ied
 25            rows(kk) = rows(kk) + omega*wksp(kk-ist+1)
           if (l == 0  .or.  icode /= 1) go to 35
                 do 30 kk = ist,ied
 30                  rows(kk+l) = rows(kk+l) + omega*wksp(kk-ist+1)
 35            if (icode == 1) go to 45
           do 40 kk = ist,ied
 40            cols(kk+l) = cols(kk+l) + omega*wksp(kk-ist+1)
 45         continue
        if (j == 0) go to 90
        l1 = -i - j
        ia2 = max (ia1,1-l1)
        ib2 = min (ib1,n-l1)
        do 85 lb = 1,mb
           k = jb(lb) - incb
           l = l1 + k
           do 50 ld = 1,md
              if (jd(ld)-incd == l) go to 55
 50            continue
           go to 85
 55            ist = max (ia2,1-l)
           ied = min (ib2,n-l)
           do 60 kk = ist,ied
 60            wksp(kk-ist+1) = c(kk-i,lc)*a(kk+l1,la)*b(kk+l1,lb)
           do 65 kk = ist,ied
 65            rows(kk) = rows(kk) + omega*wksp(kk-ist+1)
           if (l == 0  .or.  icode /= 1) go to 75
                 do 70 kk = ist,ied
 70                  rows(kk+l) = rows(kk+l) + omega*wksp(kk-ist+1)
 75            if (icode == 1) go to 85
           do 80 kk = ist,ied
 80            cols(kk+l) = cols(kk+l) + omega*wksp(kk-ist+1)
 85         continue
 90      continue
 95   continue
  return
end
subroutine tsumn (nn,np,nq,lda,ldb,ldc,ldj,ma,mb,mc,md,incb,incc,incd, &
  ja,jb,jc,jd,a,b,c,rows,omega)
!
!*******************************************************************************
!
!! TSUMN computes the row sum of C*A*B restricted to the sparsity pattern of JD.
!
!
!      c is np x nn        b is nn x nq
!      a is nn x nn        d is np x nq
!
!  definition of parameters --
!
!         nn,np,nq      orders of arrays
!         lda,ldb,ldc   row dimensions of arrays a,b,c
!         ldj           row dimension of ja,jb,jc,jd vectors
!         ma,mb,mc,md   columns (diagonals) in arrays a,b,c,d
!         incb,incc,    offsets for diagonal numbers of b,c,d arrays
!           incd
!         ja,jb,jc,jd   diagonal index arrays for a,b,c,d
!         a,b,c         arrays of dimension n x (ma,mb,md)
!         rows          row sum of d = c*a*b upon output
!         omega         relaxation factor between 0 and 1
!
!  
!
  integer   ja(ldj,1), jb(ldj,1), jc(ldj,1), jd(ldj,1)
  dimension a(lda,1), b(ldb,1), c(ldc,1), rows(1)
!
  n = nn
  do 40 lc = 1,mc
     i = jc(1,lc) - incc
     ia1 = max (1,1-i)
     ib1 = min (np,n-i)
     do 35 la = 1,ma
        j = ja(1,la)
        l1 = i + j
        ia2 = max (ia1,1-l1)
        ib2 = min (ib1,n-l1)
        do 30 lb = 1,mb
           k = jb(1,lb) - incb
           l = l1 + k
           do 15 ld = 1,md
              if (jd(1,ld)-incd == l) go to 20
 15            continue
           go to 30
 20            ist = max (ia2,1-l)
           ied = min (ib2,nq-l)
           do 25 m = ist,ied
 25            rows(m) = rows(m) + omega*c(m,lc)*a(m+i,la)*b(m+l1,lb)
 30         continue
 35      continue
 40   continue
  return
end
subroutine unpmdg (ndim,nn,maxnz,jcoef,coef,ncol,nc,p,ip,maxd,maxnew,jcnew, &
  wksp,iwksp,isym)
!
!*******************************************************************************
!
!! UNPMDG reverses the permutation done by PMDIAG.  
!
!
!  it permutes the matrix according to index vector ip.
!     the permuted matrix is stored in diagonal format.
!
!  Parameters:
!
!        ndim      row dimension of coef and jcoef arrays
!                   in defining routine
!        n         order of system (active row size of coef and jcoef)
!        maxnz     active column size of coef and jcoef
!        jcoef     integer array of column numbers
!        coef      real array of coefficients
!        ncolor    number of colors in the permutation (= ncol)
!        nc        integer vector of length ncolor giving the
!                   number of nodes for each color
!        p         permutation vector
!        ip        inverse permuation vector
!        maxd      active column size of permuted matrix
!        jcnew     integer array of size ncolor*max(maxnew(i))
!                   giving the diagonal numbers for each color
!        wksp      real workspace of length n
!        iwksp     integer workspace of length 2*n
!        isym      symmetric storage switch
!                   = 0    symmetric storage
!                   = 1    nonsymmetric storage
!
!  
!
  integer jcoef(2), nc(1), p(1), jcnew(ncol,1), maxnew(1), iwksp(1), ip(1)
  dimension coef(ndim,1), wksp(1)
!
!
  n = nn
  ncolor = ncol
!
!  set up pointer vector.
!
  do 10 j = 1,maxnz
     jcol = jcoef(j)
     iwksp(n+jcol) = j
 10   continue
!
!  permute rows of matrix first.
!
  do 15 j = 1,maxd
     do 12 i = 1,n
 12      wksp(i) = coef(i,j)
     call vscatr (n,wksp,ip,coef(1,j))
 15   continue
!
!  rearrange rows.
!
  ist = 1
  do 35 k = 1,ncolor
     ncc = nc(k)
     ied = ist + ncc - 1
     lim = maxnew(k)
     do 30 i = ist,ied
        iip = ip(i)
        do 20 j = 2,maxd
           wksp(j) = coef(iip,j)
           coef(iip,j) = 0.0
 20         continue
        do 25 j = 2,lim
           if (wksp(j) == 0.0) go to 25
           jcol = ip(i + jcnew(k,j)) - iip
           l = iwksp(n+jcol)
           coef(iip,l) = wksp(j)
 25         continue
 30      continue
     ist = ist + ncc
 35   continue
!
!  zero out lower triangular matrix if symmeteric storage used.
!
  if (isym /= 0) return
  maxold = (maxnz + 1)/2
  mp1 = maxold + 1
  do 45 j = mp1,maxnz
     do 40 i = 1,n
 40      coef(i,j) = 0.0
     jcoef(j) = 0
 45   continue
  maxnz = maxold
  return
end
subroutine uscal1 (nn,ndim,maxnzz,jcoef,coef,rhs,u,ubar,diag,work,iflag)
!
!*******************************************************************************
!
!! USCAL1 reverses the scaling done in routine SCAL1.  
!
!
!  diag must contain upon input the output from scal1.
!     (Purdue data structure)
!
!  Parameters:
!
!         n       dimension of matrix
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         ubar    exact solution (optional)
!         diag    vector of the same name from scal1 routine
!         work    work array of length n (volatile)
!         iflag   flag for ubar
!                  = 0  do not unscale ubar
!                  = 1  unscale ubar
!
!  
!
  integer   jcoef(ndim,1)
  dimension coef(ndim,1), rhs(1), u(1), diag(1), work(1), ubar(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
  maxnz = maxnzz
!
!  unscale u and rhs arrays.
!
  do 10 i = 1,n
 10   u(i) = diag(i)*u(i)
  if (iflag == 0) go to 20
  do 15 i = 1,n
 15   ubar(i) = diag(i)*ubar(i)
 20   do 25 i = 1,n
 25   diag(i) = 1.0/diag(i)
  do 30 i = 1,n
 30   rhs(i) = diag(i)*rhs(i)
!
!  unscale matrix.
!
  if (keygs == 2) go to 45
!
!  using gathers.
!
  do 40 j = 1,maxnz
     call vgathr (n,diag,jcoef(1,j),work)
     do 35 i = 1,n
 35      coef(i,j) = diag(i)*coef(i,j)*work(i)
 40   continue
  return
!
!  not using gathers.
!
 45   do 55 j = 1,maxnz
     do 50 i = 1,n
 50      coef(i,j) = diag(i)*coef(i,j)*diag(jcoef(i,j))
 55   continue
  return
end
subroutine uscal2 (nn,ndim,maxnz,jcoef,coef,rhs,u,ubar,diag,iflag)
!
!*******************************************************************************
!
!! USCAL2 reverses the scaling done in routine SCAL2.  
!
!
!  diag must contain upon input the output from scal2.
!     (diagonal data structure)
!
!  Parameters:
!
!         n       dimension of matrix
!         ndim    row dimension of coef array in defining routine
!         maxnz   number of columns in coef array
!         jcoef   integer matrix representation array
!         coef    matrix representation array
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         ubar    exact solution (optional)
!         diag    vector of the same name from scal2 routine
!         iflag   flag for ubar
!                  = 0  do not unscale ubar
!                  = 1  unscale ubar
!
!  
!
  integer   jcoef(2)
  dimension coef(ndim,1), rhs(1), u(1), diag(1), ubar(1)
!
!
  n = nn
!
!  unscale u and rhs arrays.
!
  do 10 i = 1,n
 10   u(i) = diag(i)*u(i)
  if (iflag == 0) go to 20
  do 15 i = 1,n
 15   ubar(i) = diag(i)*ubar(i)
 20   do 25 i = 1,n
 25   diag(i) = 1.0/diag(i)
  do 30 i = 1,n
 30   rhs(i) = diag(i)*rhs(i)
!
!  unscale matrix.
!
  do 50 j = 1,maxnz
     ind = jcoef(j)
     len = n - iabs(ind)
     if (ind < 0) go to 40
     do 35 i = 1,len
 35      coef(i,j) = diag(i)*coef(i,j)*diag(i+ind)
     go to 50
 40      do 45 i = 1,len
 45      coef(i-ind,j) = diag(i)*coef(i-ind,j)*diag(i-ind)
 50   continue
  return
end
subroutine uscal3 (nn,nz,ia,ja,a,rhs,u,ubar,diag,work,iflag)
!
!*******************************************************************************
!
!! USCAL3 reverses the scaling done in SCAL3.  
!
!
!  diag must contain upon input the output from scal3.
!     (sparse data structure)
!
!  Parameters:
!
!         n       dimension of matrix
!         nz      length of the vectors a, ia, and ja
!         a       vector of matrix coefficients
!         ia      vector of i values
!         ja      vector of j values
!         rhs     right hand side of matrix problem
!         u       latest estimate of solution
!         ubar    exact solution (optional)
!         diag    vector of the same name from scal3 routine
!         work    work array of length n (volatile)
!         iflag   flag for ubar
!                  = 0  do not unscale ubar
!                  = 1  unscale ubar
!
!  
!
  integer   ia(1), ja(1)
  dimension a(1), rhs(1), u(1), diag(1), work(1),ubar(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
!
!  unscale u and rhs arrays.
!
  do 10 i = 1,n
 10   u(i) = diag(i)*u(i)
  if (iflag == 0) go to 20
  do 15 i = 1,n
 15   ubar(i) = diag(i)*ubar(i)
 20   do 25 i = 1,n
 25   diag(i) = 1.0/diag(i)
  do 30 i = 1,n
 30   rhs(i) = diag(i)*rhs(i)
!
!  unscale matrix.
!
  if (keygs == 2) go to 50
!
!  using gathers.
!
  ist = 1
 35   ied = min (ist-1+n,nz)
  if (ied < ist) return
     len = ied - ist + 1
     call vgathr (len,diag,ia(ist),work)
     do 40 i = ist,ied
 40      a(i) = a(i)*work(i-ist+1)
     call vgathr (len,diag,ja(ist),work)
     do 45 i = ist,ied
 45      a(i) = a(i)*work(i-ist+1)
  ist = ied + 1
  go to 35
!
!  not using gathers.
!
 50   do 55 i = 1,nz
 55   a(i) = a(i)*diag(ia(i))*diag(ja(i))
  return
end
subroutine uslqw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! USLQW runs the USYMLQ algorithm.  
!
!
!  see: m. a. saunders, h. d. simon
! and e. l. yip, "two conjugate-gradient-type methods for sparse
! unsymmetric linear equations, report eta-tr-18, boeing computer
! services, seattle, washington, 1984, to appear in siam journal on
! numerical analysis.
!
! note -- this routine is still not quite optimal.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2),wfac(1), jwfac(1)
  integer vect1, vect2, os
  logical uneed
  external suba, subat, subql, subqlt, subqr, subqrt
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! next, the indexing functions.
!
  indv1(i) = vect1 + mod(i,nv)*n
  indv2(i) = vect2 + mod(i,nv)*n
  indbe1(i) = ibeta1 + mod(i,os)
  indbe2(i) = ibeta2 + mod(i,os)
  indc(i) = icos + mod(i,os)
  inds(i) = isin + mod(i,os)
  indu(i) = iu + mod(i,os+1)
  indw(i) = iw + n*mod(i,os)
!
! preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel = 12
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 996
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  if (iqr) go to 995
  if (level >= 2) write (nout,496)
496   format (' usymlq')
!
! initialize the stopping test.
!
  call inithv (0)
  zdhav = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
!  associated integer variables.
!
  os = 2
  iv =  1
  nv = os
  iw =  1
  vect1 = iw + iv*n*os
  vect2 = vect1 + iv*n*nv
  ibeta1 = vect2 + iv*n*nv
  ibeta2 = ibeta1 + os
  icos = ibeta2 + os
  isin = icos + os
  iu = isin + os
  iv1 = iu + os + 1
  iv2 = iv1 + n
  nwusd = max(nwusd,iv2-1+n)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
  in = 0
  is = 0
  uneed = rcalp .or. zcalp .or. ztcalp .or. udhav .or. ntest == 6 .or. &
    level >= 3
!
!  Begin iteration loop.
!
! perform first-iterate calculations.
!
 10   if (in /= 0) go to 100
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  vnorm1 = sqrt(vdot(n,wk(iv2),wk(iv2)))
  vnorm2 = vnorm1
  if (abs(vnorm1) < srelpr) go to 997
  gamma1 = 1.0/vnorm1
  gamma2 = 1.0/vnorm2
  call vtriad (n,wk(indv1(0)),xxx,gamma1,wk(iv2),2)
  call vcopy (n,wk(indv1(0)),wk(indv2(0)))
  zdot = vnorm1**2
  ucnp1= 0.0
!
! determine whether or not to stop --
!
 100  call inithv (1)
  zdhav = .true.
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
!  compute q(n+1), etc -- the direction vectors
!
  call suba (coef,jcoef,wfac,jwfac,n,wk(indv1(in)),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  an = vdot (n,wk(indv2(in)),wk(iv2))
  if (in /= 0) go to  110
  call vtriad (n,wk(indv2(in+1)),wk(iv2),-an,wk(indv2(in)), 1)
  wk(indbe2(in)) = -an
  go to 111
 110  call vtriad (n,wk(indv2(in+1)),xxx,-vnorm1,wk(indv2(in-1)),2)
  call vtriad (n,wk(indv2(in+1)),wk(indv2(in+1)),1.0,wk(iv2),1)
  call vtriad (n,wk(indv2(in+1)),wk(indv2(in+1)),-an,wk(indv2(in)),1)
  wk(indbe2(in)) = -an
  wk(indbe2(in-1)) = -vnorm1
 111  vn2old = vnorm2
  vnorm2 = sqrt(vdot (n,wk(indv2(in+1)),wk(indv2(in+1))))
  if (abs(vnorm2) < srelpr) go to 997
  gamma2 = 1.0/vnorm2
  call vtriad (n,wk(indv2(in+1)),xxx,gamma2,wk(indv2(in+1)),2)
!
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(indv2(in)),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  if (in /= 0) go to 810
  call vtriad (n,wk(indv1(in+1)),wk(iv2),-an,wk(indv1(in)),1)
  wk(indbe1(in)) = -an
  go to 811
 810  call vtriad (n,wk(indv1(in+1)),xxx,-vn2old,wk(indv1(in-1)),2)
  call vtriad (n,wk(indv1(in+1)),wk(indv1(in+1)),1.0,wk(iv2),1)
  call vtriad (n,wk(indv1(in+1)),wk(indv1(in+1)),-an,wk(indv1(in)),1)
  wk(indbe1(in)) = -an
  wk(indbe1(in-1)) = -vn2old
 811  vn1old= vnorm1
  vnorm1 = sqrt(vdot (n,wk(indv1(in+1)),wk(indv1(in+1))))
  if (abs(vnorm1) < srelpr) go to 997
  gamma1 = 1.0/vnorm1
  call vtriad (n,wk(indv1(in+1)),xxx,gamma1,wk(indv1(in+1)),2)
!
!  now update the factorization
  ucnbar = ucnp1
  ibgn = max(0,in+1-os)
  do 1 i = ibgn,in
 1    wk(indu(i+1)) = -wk(indbe2(i))
  if (ibgn > 0) wk(indu(ibgn))= 0.0
  call qrupd (in+1,os+1,os,wk(icos),wk(isin),ucnbar,ucn,wk(iu),vn2old,ier)
  if (ier /= 0) go to 998
  ucnp1 = wk(indu(in+1))
!
!  update the old w vector.
!
  if (in /= 0)call vtriad (n,wk(indw(in-1)),xxx,ucnbar/ucn,wk(indw(in-1)),2)
!
!  now generate the new w vector.
!
  if (abs(ucnp1) < srelpr) go to 998
  call vcopy (n,wk(indv1(in)),wk(iv1))
  ibgn = max(1,in-os+1)
  iend = in
  if (iend < ibgn) go to 200
  do 201 i = ibgn,iend
 201  call vtriad (n,wk(iv1),wk(iv1),-wk(indu(i)),wk(indw(i-1)),1)
 200  continue
  call vtriad (n,wk(indw(in)),xxx,1.0/ucnp1,wk(iv1),2)
  if (in /= 0) go to 205
!
!  update iterate u(0).
  zold= 0.0
  zbar = vn1old
  if (uneed) call vtriad (n,u,u,zbar,wk(indw(0)), 1)
  go to 210
!
!  update subsequent iterates u(n).
!
 205  zold = wk(indc(in))*zbar
  zbold = zbar
  zbar =-wk(inds(in))*zbar
  factor = zold
  if (uneed) factor = factor - zbold*ucn/ucnbar
  call vtriad (n,u,u,factor,wk(indw(in-1)), 1)
  if (uneed) call vtriad (n,u,u,zbar,wk(indw(in)), 1)
 210  continue
  zdot = (zbar/ucnp1*vnorm1)**2
!
! proceed to next iteration
!
  in = in +  1
  is = is + 1
  go to  10
!
!  Finish up.
!
 900  if (.not. uneed) call vtriad (n,u,u,zbar,wk(indw(in-1)),1)
  if (halt) go to 715
  ier = 1
  call ershow (ier,'uslqw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' usymlq converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
! unimplemented option
 995  ier = -16
  call ershow (ier,'uslqw')
  return
 996  call ershow (ier,'uslqw')
  go to 735
!
 997  ier = -13
  call ershow (ier,'uslqw')
  go to 725
!
 998  ier = -14
  call ershow (ier,'uslqw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'uslqw')
  go to 735
!
end
subroutine usqrw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef,wfac, &
  jwfac,n,u,ubar,rhs,wk,nw,iparm,rparm,ier)
!
!*******************************************************************************
!
!! USQRW runs the USYMQR algorithm.
!
!
!  same reference as usymlq algorithm.
!
  dimension u(1), ubar(1), rhs(1), wk(1), coef(1), jcoef(2),wfac(1), jwfac(1)
  integer vect1, vect2, os
  external suba, subat, subql, subqlt, subqr, subqrt
  dimension iparm(30), rparm(30)
  logical iql, iqr
!
!
!
  common / itcom1 / in, itmax, level, nout, ns1, ns2, ns3, iplr, iqlr, ntest, &
         is, iacel, idgts, nbl1d, nbl2d
  logical           halt, maxadp, minadp, maxadd, minadd
  common / itcom2 / halt, maxadp, minadp, maxadd, minadd
  common / itcom3 / alpha, beta, zeta, emax, emin, pap, alphao, gamma, &
         sigma, rr, rho, dkq, dkm1, ff, rqmin, rqmax, stptst, udnm, ubarnm, &
         bnorm, bnorm1
  common / itcom4 / keygs, srelpr, keyzer
  common / itcom9 / rhave, zhave, zthave, rcalp, zcalp, ztcalp, udhav, rdhav, &
         rzhav, rzthav, zdhav, zzthav, ztdhav, rdot, rzdot, rztdot, zdot, &
         zztdot, ztdot
  logical rhave, zhave, zthave, rcalp, zcalp, ztcalp
  logical udhav, rdhav, rzhav, rzthav, zdhav, zzthav, ztdhav
!
!
!
! next, the indexing functions.
!
  indv1(i)= vect1 + mod(i,nv)*n
  indv2(i) = vect2 + mod(i,nv)*n
  indbe1(i)= ibeta1 + mod(i,os)
  indbe2(i) = ibeta2 + mod(i,os)
  indc(i) = icos + mod(i,os+ 1)
  inds(i) = isin + mod(i,os+ 1)
  indu(i) = iu + mod(i,os+2)
  indw(i) = iw + n*mod(i,os)
!
! preliminary calculations.
!
  nwusd = 0
  ier = 0
  iacel =  13
  t1 = timer (dummy)
  call echall (n,iparm,rparm,1,2,ier)
  if (ier < 0) go to 996
  iql = iqlr == 1 .or. iqlr == 3
  iqr = iqlr == 2 .or. iqlr == 3
  if (iqr) go to 995
  if (level >= 2) write (nout,496)
496   format (' usymqr')
!
! initialize the stopping test.
!
  call inithv (0)
  zdhav = .true.
  nwpstp =  nw
  call pstop (0,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk,nwpstp,ier)
  nwusd = max(nwusd,nwpstp)
  if (ier < 0) go to 730
!
!  associated integer variables.
!
  os = 2
  iv =  1
  nv = os
  iw =  1
  vect1 = iw + iv*n*os
  vect2 = vect1 + iv*n*nv
  ibeta1 = vect2 + iv*n*nv
  ibeta2 = ibeta1 + os
  icos = ibeta2 + os
  isin = icos + os+ 1
  iu = isin + os+ 1
  iv1 = iu + os+2
  iv2 = iv1 + n
  nwusd = max(nwusd,iv2- 1+n)
!
! check the memory usage --
!
  if (nwusd > nw) go to 999
!
!
! now, perform first-iterate calculations
  in = 0
  is = 0
  call suba (coef,jcoef,wfac,jwfac,n,u,wk(iv1))
  call vexopy (n,wk(iv1),rhs,wk(iv1),2)
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  vnorm1 = sqrt(vdot (n,wk(iv2),wk(iv2)))
  vnorm2 = vnorm1
  if (abs(vnorm1) < srelpr) go to 997
  gamma1 = 1.0/vnorm1
  gamma2 = 1.0/vnorm2
  call vtriad (n,wk(indv1(0)),xxx,gamma1,wk(iv2),2)
  call vcopy (n,wk(indv1(0)),wk(indv2(0)))
  zdot = vnorm1**2
  znext = vnorm1
!
!  Begin iteration loop.
!
! determine whether or not to stop --
!
 10   call inithv (1)
  zdhav = .true.
  nwpstp = nw - (iv1-1)
  call pstop (1,suba,subql,subqr,coef,jcoef,wfac,jwfac,n,u,ubar,rhs, &
    xxx,xxx,xxx,wk(iv1),nwpstp,ier)
  nwusd = max(nwusd,nwpstp+iv1-1)
  if (level >= 2) call iterm (n,u)
  if (halt .or. in >= itmax .or. ier < 0) go to 900
!
!
!  compute q(n+1), etc -- the direction vectors
  call suba (coef,jcoef,wfac,jwfac,n,wk(indv1(in)),wk(iv1))
  call subql (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  an = vdot (n,wk(indv2(in)),wk(iv2))
  if (in /= 0) go to 110
  call vtriad (n,wk(indv2(in+1)),wk(iv2),-an,wk(indv2(in)),1)
  wk(indbe2(in)) = -an
  go to 111
 110  call vtriad (n,wk(indv2(in+1)),xxx,-vnorm1,wk(indv2(in-1)),2)
  call vtriad (n,wk(indv2(in+1)),wk(indv2(in+1)),1.0,wk(iv2),1)
  call vtriad (n,wk(indv2(in+1)),wk(indv2(in+1)),-an,wk(indv2(in)),1)
  wk(indbe2(in)) = -an
  wk(indbe2(in-1)) = -vnorm1
 111  vn2old = vnorm2
  vnorm2 = sqrt(vdot (n,wk(indv2(in+1)),wk(indv2(in+1))))
  if (abs(vnorm2) < srelpr) go to 997
  gamma2 = 1.0/vnorm2
  call vtriad (n,wk(indv2(in+1)),xxx,gamma2,wk(indv2(in+1)),2)
!
  call subqlt (coef,jcoef,wfac,jwfac,n,wk(indv2(in)),wk(iv1))
  call subat (coef,jcoef,wfac,jwfac,n,wk(iv1),wk(iv2))
  if (in /= 0) go to 810
  call vtriad (n,wk(indv1(in+1)),wk(iv2),-an,wk(indv1(in)),1)
  wk(indbe1(in)) = -an
  go to 811
 810  call vtriad (n,wk(indv1(in+1)),xxx,-vn2old,wk(indv1(in-1)),2)
  call vtriad (n,wk(indv1(in+1)),wk(indv1(in+1)),1.0,wk(iv2),1)
  call vtriad (n,wk(indv1(in+1)),wk(indv1(in+1)),-an,wk(indv1(in)),1)
  wk(indbe1(in)) = -an
  wk(indbe1(in-1)) = -vn2old
 811  vnorm1 = sqrt(vdot (n,wk(indv1(in+1)),wk(indv1(in+1))))
  if (abs(vnorm1) < srelpr) go to 997
  gamma1 = 1.0/vnorm1
  call vtriad (n,wk(indv1(in+1)),xxx,gamma1,wk(indv1(in+1)),2)
!
!  now update the factorization
  ibgn = max(0,in+1-os)
  do 1 i = ibgn,in
 1    wk(indu(i+1)) = -wk(indbe2(i))
  if (ibgn > 0) wk(indu(ibgn))= 0.0
  wk(indu(in+2)) = vnorm2
  call qrupd (in+2,os+2,os+1,wk(icos),wk(isin),wk(indu(in+1)),x,wk(iu),vnorm2, &
    ier)
  if (ier < 0) go to 998
!
!  now generate the new w vector.
  uc = wk(indu(in+1))
  if (abs(uc) < srelpr) go to 998
  call vcopy (n,wk(indv1(in)),wk(iv1))
  ibgn = max(1,in-os+1)
  iend = in
  if (iend < ibgn) go to 200
  do 201 i = ibgn,iend
 201  call vtriad (n,wk(iv1),wk(iv1),-wk(indu(i)),wk(indw(i-1)),1)
 200  continue
  call vtriad (n,wk(indw(in)),xxx,1.0/uc,wk(iv1),2)
!
!  update iterates u(n).
  z = wk(indc(in+1))*znext
  znext = -wk(inds(in+1))*znext
  call vtriad (n,u,u,z,wk(indw(in)),1)
  zdot = znext**2
!
! proceed to next iteration
!
  in = in + 1
  is = is + 1
  go to 10
!
!  Finish up.
!
 900  if (halt) go to 715
  ier = 1
  call ershow (ier,'usqrw')
  zeta = stptst
  go to 725
 715  continue
  if (level >= 1) write (nout,720) in
 720  format (/' usymqr converged in ',i5,' iterations.')
!
 725  continue
  if (idgts < 0) go to 730
  call perror2 (suba,coef,jcoef,wfac,jwfac,n,u,rhs,wk,digit1,digit2,idgts)
 730  t2 = timer (dummy)
  timit = t2 - t1
  iparm(2) = in
  rparm(1) = zeta
  rparm(2) = emax
  rparm(3) = emin
  rparm(6) = timit
  rparm(7) = digit1
  rparm(8) = digit2
 735  continue
  if (level >= 3) call echall (n,iparm,rparm,2,2,ier)
  nw = nwusd
  return
!
! error returns
!
 995  ier = -16
  call ershow (ier,'usqrw')
  return
!
 996  call ershow (ier,'usqrw')
  go to 735
!
 997  ier = -13
  call ershow (ier,'usqrw')
  go to 725
!
 998  ier = -14
  call ershow (ier,'usqrw')
  go to 725
!
 999  ier = -2
  call ershow (ier,'usqrw')
  go to 735
!
end
subroutine usymlq (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef, &
  n,u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! USYMLQ is the user interface to the USYMLQ algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2),wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call uslqw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef, &
    wksp,iwksp,n,u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
subroutine usymqr (suba,subat,subql,subqlt,subqr,subqrt,subadp,coef,jcoef,n, &
  u,ubar,rhs,wksp,iwksp,iparm,rparm,ier)
!
!*******************************************************************************
!
!! USYMQR is the user interface to the USYMQR algorithm.
!
  dimension u(1), ubar(1), rhs(1), coef(1), jcoef(2),wksp(1), iwksp(1)
  dimension iparm(30), rparm(30)
  external suba,  subql,  subqr
  external subat, subqlt, subqrt
  external subadp
!
!  data common blocks
!
  common / cwkcon / lenr, irpnt, irmax, leni, iipnt, iimax
!
  nw = lenr - irpnt + 1
  call usqrw (suba,subat,subql,subqlt,subqr,subqrt,coef,jcoef, &
    wksp,iwksp,n,u,ubar,rhs,wksp(irpnt),nw,iparm,rparm,ier)
  irmax = max (irmax,irpnt-1+nw)
  iimax = max (iimax,iipnt-1)
  return
end
function vadd (n,v)
!
!*******************************************************************************
!
!! VADD adds the elements of a vector.
!
  dimension v(1)
!
  sum = 0e0
  do 1 i=1,n
 1    sum = sum + v(i)
  vadd = sum
  return
end
subroutine vaddd (lda,ldja,nn,m,mdiagg,a,ja,y,x,jofff)
!
!*******************************************************************************
!
!! VADDD computes  y = y + A*x.  (diagonal storage)
!
!
!  Parameters:
!
!         lda      leading dimension of a array
!         ldja     leading dimension of ja array
!         n        active row size of matrix
!         m        active column size of matrix
!         mdiag    number of diagonals in a
!         a        array of matrix diagonals
!         ja       array of matrix diagonal numbers
!         y,x      vectors of length n
!         joff     offset for diagonal numbers
!
!  
!
  dimension a(lda,3), x(1), y(1)
  integer   ja(ldja,3)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
  mdiag = mdiagg
  joff = jofff
  if (mdiag < 1) return
  if (keyzer == 1) go to 20
  do 15 j = 1,mdiag
     ind = ja(1,j) - joff
     ist = max (1,1-ind)
     ied = min (n,m-ind)
     do 10 i = ist,ied
 10      y(i) = y(i) + a(i,j)*x(i+ind)
 15   continue
  return
!
!  unrolled version (requires memory to be zeroed out).
!
 20   l = mod (mdiag,4)
  if (l == 0) go to 60
!
!  initial short computations
!
  go to (25,35,45), l
 25   do 30 i = 1,n
 30   y(i) = y(i) + a(i,1)*x(i+ja(1,1)-joff)
  go to 55
 35   do 40 i = 1,n
 40   y(i) = y(i) + a(i,1)*x(i+ja(1,1)-joff) + a(i,2)*x(i+ja(1,2)-joff)
  go to 55
 45   do 50 i = 1,n
 50   y(i) = y(i) + a(i,1)*x(i+ja(1,1)-joff) + a(i,2)*x(i+ja(1,2)-joff) &
        + a(i,3)*x(i+ja(1,3)-joff)
 55   if (mdiag <= 4) return
!
!  loop unrolling to a level of 4.
!
 60   lp1 = l + 1
  do 70 j = lp1,mdiag,4
     do 65 i = 1,n
 65      y(i) = y(i) + a(i,j  )*x(i+ja(1,j  )-joff) &
                     + a(i,j+1)*x(i+ja(1,j+1)-joff) &
                     + a(i,j+2)*x(i+ja(1,j+2)-joff) &
                     + a(i,j+3)*x(i+ja(1,j+3)-joff)
 70   continue
  return
end
subroutine vadddt (lda,ldja,nn,m,mdiagg,a,ja,y,x,jofff)
!
!*******************************************************************************
!
!! VADDDT computes  y = y + (A**t)*x.  (diagonal storage)
!
!
!  Parameters:
!
!         lda      leading dimension of a array
!         ldja     leading dimension of ja array
!         n        active row size of matrix
!         m        active column size of matrix
!         mdiag    number of diagonals in a
!         a        array of matrix diagonals
!         ja       array of matrix diagonal numbers
!         y,x      vectors of length n
!         joff     offset for diagonal numbers
!
!  
!
  dimension a(lda,3), x(1), y(1)
  integer   ja(ldja,3)
!
  n = nn
  mdiag = mdiagg
  joff = jofff
  if (mdiag < 1) return
  do 15 j = 1,mdiag
     ind = ja(1,j) - joff
     ist = max (1,1-ind)
     ied = min (n,m-ind)
     do 10 i = ist,ied
 10      y(i+ind) = y(i+ind) + a(i,j)*x(i)
 15   continue
  return
end
subroutine vaddp (ndimr,ndimi,nn,mm,a,ja,y,x,wksp)
!
!*******************************************************************************
!
!! VADDP does  y = y + A*x  (Purdue format)
!
!
!  Parameters:
!
!       ndimr     row dimension of a array
!       ndimi     row dimension of ja array
!       n         order of system
!       m         number of columns in a and ja arrays
!       a         real array of active size n by m
!       ja        integer array of active size n by m
!       y         accumulation vector
!       x         right-hand-side vector
!       wksp      workspace vector of length n  (keygs = 1 only)
!
!  
!
  dimension a(ndimr,3), ja(ndimi,3), x(1), y(1), wksp(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
  m = mm
  if (m <= 0) return
  if (keygs == 1) go to 100
!
!  implicit gathers.
!
  l = mod (m,4)
  if (l == 0) go to 45
!
!  initial short computations
!
  go to (10,20,30), l
 10   do 15 i = 1,n
 15   y(i) = y(i) + a(i,1)*x(ja(i,1))
  go to 40
 20   do 25 i = 1,n
 25   y(i) = y(i) + a(i,1)*x(ja(i,1)) + a(i,2)*x(ja(i,2))
  go to 40
 30   do 35 i = 1,n
 35   y(i) = y(i) + a(i,1)*x(ja(i,1)) + a(i,2)*x(ja(i,2))+ a(i,3)*x(ja(i,3))
 40   if (m <= 4) return
!
!  loop unrolling to a level of 4.
!
 45   lp1 = l + 1
  do 55 j = lp1,m,4
     do 50 i = 1,n
 50      y(i) = y(i) + a(i,j)*x(ja(i,j)) + a(i,j+1)*x(ja(i,j+1)) &
                     + a(i,j+2)*x(ja(i,j+2)) + a(i,j+3)*x(ja(i,j+3))
 55   continue
  return
!
!  explicit gathers.
!
 100  do 110 j = 1,m
     call vgathr (n,x,ja(1,j),wksp)
     do 105 i = 1,n
 105     y(i) = y(i) + a(i,j)*wksp(i)
 110  continue
  return
end
subroutine vaddpt (ndimr,ndimi,n,m,a,ja,y,x,wksp)
!
!*******************************************************************************
!
!! VADDPT does  y = y + (A**t)*x  (Purdue format)
!
!
!  Parameters:
!
!       ndimr     row dimension of a array
!       ndimi     row dimension of ja array
!       n         order of system
!       m         number of columns in a and ja arrays
!       a         real array of active size n by m
!       ja        integer array of active size n by m
!       y         accumulation vector
!       x         right-hand-side vector
!       wksp      workspace vector of length n
!
!  
!
  dimension a(ndimr,3), ja(ndimi,3), x(1), y(1), wksp(1)
!
  if (m <= 0) return
!
  do 20 j = 1,m
     do 15 i = 1,n
        y(ja(i,j)) = y(ja(i,j)) + a(i,j)*x(i)
 15      continue
 20   continue
  return
end
subroutine vadds (mm,np,ia,ja,a,y,x,wksp)
!
!*******************************************************************************
!
!! VADDS does  y = y + A*x  (sparse format)
!
!
!  Parameters:
!
!       m         number of partitions
!       np        partition pointers
!       ia        vector of i values
!       ja        vector of j values
!       a         vector of coefficients
!       y         accumulation vector
!       x         right-hand-side vector
!       wksp      workspace vector of length 2*n  (keygs = 1 only)
!
!  
!
  dimension np(1), a(1), ia(1), ja(1), x(1), y(1), wksp(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  m = mm
  if (m <= 0) return
  if (keygs == 1) go to 20
!
!  implicit gathers.
!
  do 15 k = 1,m
     ist = np(k)
     ied = np(k+1) - 1
!dir$ ivdep
     do 10 i = ist,ied
 10      y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
 15   continue
  return
!
!  explicit gathers.
!
 20   do 30 k = 1,m
     ist = np(k)
     nel = np(k+1) - ist
     call vgathr (nel,x,ja(ist),wksp)
     call vgathr (nel,y,ia(ist),wksp(nel+1))
     do 25 i = 1,nel
 25      wksp(i) = wksp(nel+i) + a(ist+i-1)*wksp(i)
     call vscatr (nel,wksp,ia(ist),y)
 30   continue
  return
end
subroutine vcopy (n,x,y)
!
!*******************************************************************************
!
!! VCOPY copies vector x to vector y.
!
!
!  Parameters:
!
!          n       length of vectors
!          x,y     vectors of length n
!
  integer n
!
  real x(n)
  real y(n)
!
  y = x

  return
end
function vdot (n,x,y)
!
!*******************************************************************************
!
!! VDOT computes the dot product of two vectors.
!
!
!  Parameters:
!
!          n       length of vectors
!          x,y     vectors of length n
!
!  
!
  dimension x(1), y(1)
!
!
  vdot = 0.0
  if (n <= 0) return
  do 10 i = 1,n
 10   vdot = vdot + x(i)*y(i)
  return
end
subroutine vemxty (nn,v,x,y)
!
!*******************************************************************************
!
!! VEMTXY computes  v = -x * y  where v, x, and y are vectors.
!
!
!  Parameters:
!
!          n       length of vectors  (= nn)
!          v,x,y   vectors of length n
!
!  
!
  dimension v(1), x(1), y(1)
!
  n = nn
  if (n <= 0) return
  do 10 i = 1,n
 10   v(i) = -x(i)*y(i)
  return
end
subroutine vexopy (nn,v,x,y,icode)
!
!*******************************************************************************
!
!! VEXOPY computes  v = x  op  y  where  op is one of the operations  + - * / .
!
!
!  Parameters:
!
!          n       length of vectors  (= nn)
!          v,x,y   vectors of length n
!          icode   key indicating operation
!            = 1      for addition
!            = 2      for subtraction
!            = 3      for multiplication
!            = 4      for division
!
!  
!
  dimension v(1), x(1), y(1)
!
  n = nn
  if (n <= 0) return
  go to (10,20,30,40), icode
!
!  compute   v = x + y
!
 10   do 15 i = 1,n
 15   v(i) = x(i) + y(i)
  return
!
!  compute   v = x - y
!
 20   do 25 i = 1,n
 25   v(i) = x(i) - y(i)
  return
!
!  compute   v = x * y
!
 30   do 35 i = 1,n
 35   v(i) = x(i)*y(i)
  return
!
!  compute   v = x / y
!
 40   do 45 i = 1,n
 45   v(i) = x(i)/y(i)
  return
end
subroutine vfill (n,v,val)
!
!*******************************************************************************
!
!! VFILL fills a vector with a constant value.
!
!
!  Parameters:
!
!          n      integer length of vector v
!          v      vector
!          val    constant that fills first n locations of v
!
!  
!
  dimension v(n)
!
  if (n <= 0) return
  do 10 i = 1,n
 10   v(i) = val
  return
end
subroutine vgathi (n,ja,ia,jb)
!
!*******************************************************************************
!
!! VGATHI gathers elements from an array.
!
!
!  The elements are gathered according to index
!  list ia and places them into consecutive locations in
!  array jb.
!
!  Parameters:
!
!          n        order of arrays ia and jb
!          ja       integer array of source elements
!          ia       integer array of length n giving desired
!                      elements of array ja
!          jb       integer target array of length n
!
!  
!
  integer   ia(1), ja(1), jb(1)
!
  if (n <= 0) return
  do 10 i = 1,n
 10   jb(i) = ja(ia(i))
!
!205  jb(1;n) = q8vgathr (ja(1;n),ia(1;n);jb(1;n))
!ray1 call gather (n,jb,ja,ia)
!
  return
end
subroutine vgathr (n,a,ia,b)
!
!*******************************************************************************
!
!! VGATHR gathers elements from an array.
!
!
!  The elements are gathered according to index
!  list ia and places them into consecutive locations in
!  array b.
!
!  Parameters:
!
!          n       order of arrays ia and b
!          a       array of source elements
!          ia      integer array of length n giving desired
!                     elements of array a
!          b       target array of length n
!
!  
!
  integer   ia(1)
  dimension a(1), b(1)
!
  if (n <= 0) return
  do 10 i = 1,n
 10   b(i) = a(ia(i))
!
!205  b(1;n) = q8vgathr (a(1;n),ia(1;n);b(1;n))
!ray1 call gather (n,b,a,ia)
!
  return
end
subroutine vicopy (n,iv1,iv2)
!
!*******************************************************************************
!
!! VICOPY copies one integer vector to another.
!
  integer iv1(1), iv2(1)
  if (n <= 0) return
  do 1 i=1,n
 1    iv2(i) = iv1(i)
  return
end
subroutine vifill (n, iv, ival)
!
!*******************************************************************************
!
!! VIFILL fills an integer vector with a value.
!
  integer iv(1)
  if (n <= 0) return
  do 1 i=1,n
 1    iv(i) = ival
  return
end
subroutine vinv ( n, v )
!
!*******************************************************************************
!
!! VINV computes v = 1/v.
!
!
!  Parameters:
!
!        n       length of vector
!        v       input/output vector of length n.
!
  integer n
!
  integer i
  real v(n)
!
  do i = 1, n
    v(i) = 1.0 / v(i)
  end do

  return
end
function vmax (n,v)
!
!*******************************************************************************
!
!! VMAX determines the maximum element of a vector.
!
!
!  Parameters:
!
!        n     length of vector
!        v     real vector of length n
!
!  
!
  dimension v(1)
!
  vmax = v(1)
  if (n <= 1) return
  do 10 i = 2,n
     if (v(i) > vmax) vmax = v(i)
 10   continue
!205  vmax = q8smax (v(1;n))
!ray  imax = ismax (n,v,1)
!ray  vmax = v(imax)
  return
end
function vmin (n,v)
!
!*******************************************************************************
!
!! VMIN determines the minimum element of a vector v.
!
!
!  Parameters:
!
!        n     length of vector
!        v     real vector of length n
!
!  
!
  dimension v(1)
!
  vmin = v(1)
  if (n <= 1) return
  do 10 i = 2,n
     if (v(i) < vmin) vmin = v(i)
 10   continue
!205  vmin = q8smin (v(1;n))
!ray  imin = ismin (n,v,1)
!ray  vmin = v(imin)
  return
end
subroutine vscati (n,ja,ia,jb)
!
!*******************************************************************************
!
!! VSCATI scatters elements from consecutive locations in an array.
!
!
!  The elements are scattered to positions in array jb according to index list ia.
!
!  Parameters:
!
!         n       order of arrays ia and ja
!         ja      integer array of source elements
!         ia      integer array of length n giving new locations
!                   in array jb.
!         jb      integer target array
!
!  
!
  integer   ia(1), ja(1), jb(1)
!
  if (n <= 0) return
  do 10 i = 1,n
 10   jb(ia(i)) = ja(i)
!
!205  jb(1;n) = q8vscatr (ja(1;n),ia(1;n);jb(1;n))
!ray1 call scatter (n,jb,ia,ja)
!
  return
end
subroutine vscatr (n,a,ia,b)
!
!*******************************************************************************
!
!! VSCATR scatters elements from consecutive locations in an array.
!
!
!  The elements are scattered to positions in array b according to index list ia.
!
!  Parameters:
!
!          n       order of arrays ia and a
!          a       array of source elements
!          ia      integer array of length n giving new locations
!                    in array b
!          b       target array
!
!  
!
  integer   ia(1)
  dimension a(1), b(1)
!
  if (n <= 0) return
  do 10 i = 1,n
 10   b(ia(i)) = a(i)
!
!205  b(1;n) = q8vscatr (a(1;n),ia(1;n);b(1;n))
!ray1 call scatter (n,b,ia,a)
!
  return
end
subroutine vsqrt (n,v,w)
!
!*******************************************************************************
!
!! VSQRT computes the square root of the entries of a vector.
!
!
  dimension v(1), w(1)
!
  do 1 i=1,n
 1    w(i) = sqrt(v(i))
  return
end
subroutine vsrta1 (nz,ia,ja,a)
!
!*******************************************************************************
!
!! VSRTA1 sorts the sparse data structure by rows and then columns.
!
!
!   usage               - call vsrta1 (nz,ia,ja,a)
!
!   arguments    ia     - on input, ia contains the i values of
!                          the array to be sorted.
!                         on output, ia contains the i values of
!                          the sorted array.
!                ja     - on input, ja contains the j values of
!                          the array to be sorted.
!                         on output, ja contains the j values of
!                          the sorted array.
!                a      - on input, a contains the coefficients
!                          of the array to be sorted.
!                         on output, a contains the coefficients
!                          of the sorted array.
!                nz     - input variable containing the number of
!                           elements in the array to be sorted.
!
!   precision/hardware  - single/all
!
!   reqd. imsl routines - none required
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp
!
!   copyright           - 1978 by imsl, inc. all rights reserved.
!
!   warranty            - imsl warrants only that imsl testing has been
!                           applied to this code.  no other warranty,
!                           expressed or implied, is applicable.
!
  integer            ia(nz), ja(nz)
  dimension          a(nz)
!                                  specifications for local variables
  integer            iu(21),il(21)
!
  logical lt, le, eq
  lt (i1,j1,i2,j2) = i1<i2 .or. (i1==i2 .and. j1<j2)
  le (i1,j1,i2,j2) = i1<i2 .or. (i1==i2 .and. j1<=j2)
  eq (i1,j1,i2,j2) = i1==i2 .and. j1==j2
!
  m = 1
  i = 1
  j = nz
  r = 0.375
  if (nz <= 0) return
   10 if (i == j) go to 55
  if (r > 0.5898437) go to 20
  r = r + 3.90625e-2
  go to 25
   20 r = r - 0.21875
   25 k = i
!                                  select a central element of the
!                                  array and save it in location t
  ij = int ( float(i) + float(j-i)*r )
  t = a(ij)
  it = ia(ij)
  jt = ja(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
  if ( le(ia(i),ja(i),it,jt) ) go to 30
  ia(ij) = ia(i)
  ia(i) = it
  it = ia(ij)
  ja(ij) = ja(i)
  ja(i) = jt
  jt = ja(ij)
  a(ij) = a(i)
  a(i) = t
  t = a(ij)
   30 l = j
!                                  if last element of array is less than
!                                  t, interchange with t
  if (.not. lt(ia(j),ja(j),it,jt) ) go to 40
  ia(ij) = ia(j)
  ia(j) = it
  it = ia(ij)
  ja(ij) = ja(j)
  ja(j) = jt
  jt = ja(ij)
  a(ij) = a(j)
  a(j) = t
  t = a(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
  if ( le(ia(i),ja(i),it,jt) ) go to 40
  ia(ij) = ia(i)
  ia(i) = it
  it = ia(ij)
  ja(ij) = ja(i)
  ja(i) = jt
  jt = ja(ij)
  a(ij) = a(i)
  a(i) = t
  t = a(ij)
  go to 40
   35 if ( eq(ia(l),ja(l),ia(k),ja(k)) ) go to 40
  itt = ia(l)
  ia(l) = ia(k)
  ia(k) = itt
  jtt = ja(l)
  ja(l) = ja(k)
  ja(k) = jtt
  tt = a(l)
  a(l) = a(k)
  a(k) = tt
!                                  find an element in the second half of
!                                  the array which is smaller than t
   40 l = l - 1
  if (.not. le (ia(l),ja(l),it,jt) ) go to 40
!                                  find an element in the first half of
!                                  the array which is greater than t
   45 k = k + 1
  if ( lt (ia(k),ja(k),it,jt) ) go to 45
!                                  interchange these elements
  if (k <= l) go to 35
!                                  save upper and lower subscripts of
!                                  the array yet to be sorted
  if (l-i <= j-k) go to 50
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  go to 60
   50 il(m) = k
  iu(m) = j
  j = l
  m = m + 1
  go to 60
!                                  begin again on another portion of
!                                  the unsorted array
   55 m = m - 1
  if (m == 0) return
  i = il(m)
  j = iu(m)
   60 if (j-i >= 11) go to 25
  if (i == 1) go to 10
  i = i - 1
   65 i = i + 1
  if (i == j) go to 55
  it = ia(i+1)
  jt = ja(i+1)
  t = a(i+1)
  if ( le (ia(i),ja(i),it,jt) ) go to 65
  k = i
   70 ia(k+1) = ia(k)
  ja(k+1) = ja(k)
  a(k+1) = a(k)
  k = k - 1
  if ( lt (it,jt,ia(k),ja(k)) ) go to 70
  ia(k+1) = it
  ja(k+1) = jt
  a(k+1) = t
  go to 65
end
subroutine vsubd (lda,ldja,nn,m,mdiagg,a,ja,y,x,jofff)
!
!*******************************************************************************
!
!! VSUBD computes  y = y - A*x.  (diagonal storage)
!
!
!  Parameters:
!
!         lda      leading dimension of a array
!         ldja     leading dimension of ja array
!         n        active row size of matrix
!         m        active column size of matrix
!         mdiag    number of diagonals in a
!         a        array of matrix diagonals
!         ja       array of matrix diagonal numbers
!         y,x      vectors of length n
!         joff     offset for diagonal numbers
!
!  
!
  dimension a(lda,3), x(1), y(1)
  integer   ja(ldja,3)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
  mdiag = mdiagg
  joff = jofff
  if (mdiag < 1) return
  if (keyzer == 1) go to 20
  do 15 j = 1,mdiag
     ind = ja(1,j) - joff
     ist = max (1,1-ind)
     ied = min (n,m-ind)
     do 10 i = ist,ied
 10      y(i) = y(i) - a(i,j)*x(i+ind)
 15   continue
  return
!
!  unrolled version (requires memory to be zeroed out).
!
 20   l = mod (mdiag,4)
  if (l == 0) go to 60
!
!  initial short computations
!
  go to (25,35,45), l
 25   do 30 i = 1,n
 30   y(i) = y(i) - a(i,1)*x(i+ja(1,1)-joff)
  go to 55
 35   do 40 i = 1,n
 40   y(i) = y(i) - a(i,1)*x(i+ja(1,1)-joff) - a(i,2)*x(i+ja(1,2)-joff)
  go to 55
 45   do 50 i = 1,n
 50   y(i) = y(i) - a(i,1)*x(i+ja(1,1)-joff) - a(i,2)*x(i+ja(1,2)-joff) &
                - a(i,3)*x(i+ja(1,3)-joff)
 55   if (mdiag <= 4) return
!
!  loop unrolling to a level of 4.
!
 60   lp1 = l + 1

  do j = lp1,mdiag,4
     do i = 1,n
       y(i) = y(i) - a(i,j  )*x(i+ja(1,j  )-joff) &
                    - a(i,j+1)*x(i+ja(1,j+1)-joff) &
                    - a(i,j+2)*x(i+ja(1,j+2)-joff) &
                    - a(i,j+3)*x(i+ja(1,j+3)-joff)
     end do
  end do

  return
end
subroutine vsubdt (lda,ldja,nn,m,mdiagg,a,ja,y,x,jofff)
!
!*******************************************************************************
!
!! VSUBDT computes  y = y - (A**t)*x.  (diagonal storage)
!
!
!  Parameters:
!
!         lda      leading dimension of a array
!         ldja     leading dimension of ja array
!         n        active row size of matrix
!         m        active column size of matrix
!         mdiag    number of diagonals in a
!         a        array of matrix diagonals
!         ja       array of matrix diagonal numbers
!         y,x      vectors of length n
!         joff     offset for diagonal numbers
!
!  
!
  dimension a(lda,3), x(1), y(1)
  integer   ja(ldja,3)
!
  n = nn
  mdiag = mdiagg
  joff = jofff
  if (mdiag < 1) return
  do 15 j = 1,mdiag
     ind = ja(1,j) - joff
     ist = max (1,1-ind)
     ied = min (n,m-ind)
     do 10 i = ist,ied
 10      y(i+ind) = y(i+ind) - a(i,j)*x(i)
 15   continue
  return
end
subroutine vsubp (ndimr,ndimi,nn,mm,a,ja,y,x,wksp)
!
!*******************************************************************************
!
!! VSUBP does  y = y - A*x  (Purdue format).
!
!
!  Parameters:
!
!       ndimr     row dimension of a array
!       ndimi     row dimension of ja array
!       n         order of system
!       m         number of columns in a and ja arrays
!       a         real array of active size n by m
!       ja        integer array of active size n by m
!       y         accumulation vector
!       x         right-hand-side vector
!       wksp      workspace vector of length n  (keygs = 1 only)
!
!  
!
  dimension a(ndimr,3), ja(ndimi,3), x(1), y(1), wksp(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  n = nn
  m = mm
  if (m <= 0) return
  if (keygs == 1) go to 100
!
!  implicit gathers.
!
  l = mod (m,4)
  if (l == 0) go to 45
!
!  initial short computations
!
  go to (10,20,30), l
 10   do 15 i = 1,n
 15   y(i) = y(i) - a(i,1)*x(ja(i,1))
  go to 40
 20   do 25 i = 1,n
 25   y(i) = y(i) - a(i,1)*x(ja(i,1)) - a(i,2)*x(ja(i,2))
  go to 40
 30   do 35 i = 1,n
 35   y(i) = y(i) - a(i,1)*x(ja(i,1)) - a(i,2)*x(ja(i,2))- a(i,3)*x(ja(i,3))
 40   if (m <= 4) return
!
!  loop unrolling to a level of 4.
!
 45   lp1 = l + 1
  do 55 j = lp1,m,4
     do 50 i = 1,n
 50      y(i) = y(i) - a(i,j)*x(ja(i,j)) - a(i,j+1)*x(ja(i,j+1)) &
                    - a(i,j+2)*x(ja(i,j+2)) - a(i,j+3)*x(ja(i,j+3))
 55   continue
  return
!
!  explicit gathers.
!
 100  do 110 j = 1,m
     call vgathr (n,x,ja(1,j),wksp)
     do 105 i = 1,n
 105     y(i) = y(i) - a(i,j)*wksp(i)
 110  continue
  return
end
subroutine vsubpt (ndimr,ndimi,n,m,a,ja,y,x,wksp)
!
!*******************************************************************************
!
!! VSUBPT does  y = y - (A**t)*x  (Purdue format).
!
!
!  Parameters:
!
!       ndimr     row dimension of a array
!       ndimi     row dimension of ja array
!       n         order of system
!       m         number of columns in a and ja arrays
!       a         real array of active size n by m
!       ja        integer array of active size n by m
!       y         accumulation vector
!       x         right-hand-side vector
!       wksp      workspace vector of length n
!
!  
!
  dimension a(ndimr,3), ja(ndimi,3), x(1), y(1), wksp(1)
!
  if (m <= 0) return
!
  do 20 j = 1,m
     do 15 i = 1,n
        y(ja(i,j)) = y(ja(i,j)) - a(i,j)*x(i)
 15      continue
 20   continue
  return
end
subroutine vsubs (mm,np,ia,ja,a,y,x,wksp)
!
!*******************************************************************************
!
!! VSUBS does  y = y - A*x  (sparse format).
!
!
!  Parameters:
!
!       m         number of partitions
!       np        partition pointers
!       ia        vector of i values
!       ja        vector of j values
!       a         vector of coefficients
!       y         accumulation vector
!       x         right-hand-side vector
!       wksp      workspace vector of length 2*n  (keygs = 1 only)
!
  dimension np(1), a(1), ia(1), ja(1), x(1), y(1), wksp(1)
!
!
!
  common / itcom4 / keygs, srelpr, keyzer
!
!
!
  m = mm
  if (m <= 0) return
  if (keygs == 1) go to 20
!
!  implicit gathers.
!
  do 15 k = 1,m
     ist = np(k)
     ied = np(k+1) - 1
!dir$ ivdep
     do 10 i = ist,ied
 10      y(ia(i)) = y(ia(i)) - a(i)*x(ja(i))
 15   continue
  return
!
!  explicit gathers.
!
 20   do 30 k = 1,m
     ist = np(k)
     nel = np(k+1) - ist
     call vgathr (nel,x,ja(ist),wksp)
     call vgathr (nel,y,ia(ist),wksp(nel+1))
     do 25 i = 1,nel
 25      wksp(i) = wksp(nel+i) - a(ist+i-1)*wksp(i)
     call vscatr (nel,wksp,ia(ist),y)
 30   continue
  return
end
subroutine vtriad (n,c,a,con,b,icode)
!
!*******************************************************************************
!
!! VTRIAD computes C = A + CON*B or C = CON*B.
!
!
!  Parameters:
!
!        n         length of vectors
!        c,a,b     vectors of length n
!        con       multiplicative constant
!        icode     switch
!                  1 means compute A + CON * B
!                  2 menas compute CON * B
!
!  
!
  dimension a(1), b(1), c(1)
!
  if (n <= 0) return
  if (icode == 2) go to 15
!
!  compute    c = a + con*b
!
  do 10 i = 1,n
 10   c(i) = a(i) + con*b(i)
  return
!
!  compute    c = con*b
!
 15   do 20 i = 1,n
 20   c(i) = con*b(i)
  return
end
subroutine zbrent (n,tri,eps,nsig,aa,bb,maxfnn,ier)
!
!*******************************************************************************
!
!! ZBRENT finds a zero of a function in a change of sign interval.
!
!
!   computer            - cdc/single
!
!   latest revision     - january 1, 1978
!
!   purpose             - zero of a function which changes sign in a
!                           given interval (brent algorithm)
!
!   usage               - call zbrent (f,eps,nsig,a,b,maxfn,ier)
!
!   arguments    tri    - a tridiagonal matrix of order n
!                eps    - first convergence criterion (input).  a root,
!                           b, is accepted if abs(f(b)) is less than or
!                           equal to eps.  eps may be set to zero.
!                nsig   - second convergence criterion (input).  a root,
!                           b, is accepted if the current approximation
!                           agrees with the true solution to nsig
!                           significant digits.
!                a,b    - on input, the user must supply two points, a
!                           and b, such that f(a) and f(b) are opposite
!                           in sign. (= aa, bb)
!                           on output, both a and b are altered.  b
!                           will contain the best approximation to the
!                           root of f. see remark 1.
!                maxfn  - on input, maxfn should contain an upper bound
!                           on the number of function evaluations
!                           required for convergence.  on output, maxfn
!                           will contain the actual number of function
!                           evaluations used. (= maxfnn)
!                ier    - error parameter. (output)
!                         terminal error
!                           ier = 3 indicates the algorithm failed to
!                             converge in maxfn evaluations.
!                           ier = 4 indicates f(a) and f(b) have the
!                             same sign.
!
!   precision/hardware  - single and double/h32
!                       - single/h36,h48,h60
!
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp
!
!   remarks  1.  let f(x) be the characteristic function of the matrix
!                tri evaluated at x. function determ evaluates f(x).
!                on exit from zbrent, when ier=0, a and b satisfy the
!                following,
!                f(a)*f(b) <= 0,
!                abs(f(b)) <= abs(f(a)), and
!                either abs(f(b)) <= eps or
!                abs(a-b) <= max(abs(b),0.1)*10.0**(-nsig).
!                the presence of 0.1 in this error criterion causes
!                leading zeroes to the right of the decimal point to be
!                counted as significant digits. scaling may be required
!                in order to accurately determine a zero of small
!                magnitude.
!            2.  zbrent is guaranteed to reach convergence within
!                k = (alog((b-a)/d)+1.0)**2 function evaluations where
!                  d=min(over x in (a,b) of
!                    max(abs(x),0.1)*10.0**(-nsig)).
!                this is an upper bound on the number of evaluations.
!                rarely does the actual number of evaluations used by
!                zbrent exceed sqrt(k). d can be computed as follows,
!                  p = min (abs(a),abs(b))
!                  p = max (0.1,p)
!                  if ((a-0.1)*(b-0.1)<0.0) p = 0.1
!                  d = p*10.0**(-nsig)
!
!   copyright           - 1977 by imsl, inc. all rights reserved.
!
!   warranty            - imsl warrants only that imsl testing has been
!                           applied to this code. no other warranty,
!                           expressed or implied, is applicable.
!
  dimension          tri(2,1)
!
!
!  local package references --
!
!          determ
!                                  first executable statement
  a = aa
  b = bb
  maxfn = maxfnn
  t = 0.1**nsig
  ic = 2
  fa = determ(n,tri,a)
  fb = determ(n,tri,b)
  s = b
!                                  test for same sign
  if (fa*fb > 0.0) go to 50
    5 c = a
  fc = fa
  d = b - c
  e = d
   10 if (abs (fc) >= abs (fb)) go to 15
  a = b
  b = c
  c = a
  fa = fb
  fb = fc
  fc = fa
   15 continue
  tol = t * amax1 (abs (b),0.1)
  rm = (c - b)/2.0
!                                  test for first convergence criteria
  if (abs (fb) <= eps) go to 40
!                                  test for second convergence criteria
  if (abs (c-b) <= tol) go to 40
!                                  check evaluation counter
  if (ic >= maxfn) go to 45
!                                  is bisection forced
  if (abs (e) < tol) go to 30
  if (abs (fa) <= abs (fb)) go to 30
  s = fb/fa
  if (a /= c) go to 20
!                                  linear interpolation
  p = (c - b)*s
  q = 1.0 - s
  go to 25
!                                  inverse quadratic interpolation
   20 q = fa/fc
  r = fb/fc
  rone = r - 1.0
  p = s*((c - b)*q*(q - r) - (b - a)*rone)
  q = (q - 1.0)*rone*(s - 1.0)
   25 if (p > 0.0) q = -q
  if (p < 0.0) p = -p
  s = e
  e = d
!                                  if abs(p/q)>=75*abs(c-b) then
!                                     force bisection
  if (p + p >= 3.0*rm*q) go to 30
!                                  if abs(p/q)>=.5*abs(s) then force
!                                     bisection. s = the value of p/q
!                                     on the step before the last one
  if (p + p >= abs (s*q)) go to 30
  d = p/q
  go to 35
!                                  bisection
   30 e = rm
  d = e
!                                  increment b
   35 a = b
  fa = fb
  temp = d
  if (abs (temp) <= tol/2.0) temp = sign (tol/2.0,rm)
  b = b + temp
  s = b
  fb = determ(n,tri,s)
  ic = ic + 1
  if (fb*fc <= 0.0) go to 10
  go to 5
!                                  convergence of b
   40 a = c
  maxfn = ic
  go to 9000
!                                  maxfn evaluations
   45 ier = 3
  a = c
  maxfn = ic
  call ershow (ier,'zbrent')
  go to 9000
!                                  terminal error - f(a) and f(b) have
!                                  the same sign
   50 ier = 4
  maxfn = ic
  call ershow (ier,'zbrent')
 9000 continue
  aa = a
  bb = b
  maxfnn = maxfn
  return
end
