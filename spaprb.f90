!  spaprb.f90  17 June 2000
!
program spaprb
!
!*******************************************************************************
!
!! SPAPRB runs the SPARSPAK tests.
!
  write ( *, * ) ' '
  write ( *, * ) 'SPAPRB'
  write ( *, * ) '  Sample problems for SPARSPAK,'
  write ( *, * ) '  the Waterloo sparse matrix package.'
  write ( *, * ) ' '

  call test01
  call test02
  call test03
  call test04
  call test05
  call test06

  write ( *, * ) ' '
  write ( *, * ) 'SPAPRB'
  write ( *, * ) '  Normal end of SPARSPAK tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests the 1WD method.
!
!
!  This example involves the following equations:
!
!  2*x1  - x10       = 0
!  2*x2  - x9  - x10 = 0
!  2*x3  - x8  - x9  = 0
!  2*x4  - x7  - x8  = 0
!  2*x5  - x6  - x7  = 0
!  2*x6  - x5        = 11
!  2*x7  - x4  - x5  = 0
!  2*x8  - x3  - x4  = 0
!  2*x9  - x2  - x3  = 0
!  2*x10 - x1  - x2  = 0
!
!  with solution (1,3,5,7,9,10,8,6,4,2).
!
  integer, parameter :: maxadj = 300
  integer, parameter :: maxblk = 10
  integer, parameter :: maxenv = 300
  integer, parameter :: maxnon = 300
  integer, parameter :: n = 10
!
  integer iadj(maxadj)
  real diag(n)
  real env(maxenv)
  integer father(n)
  integer first(n)
  integer i
  integer iband
  integer ienv
  integer ierror
  integer invprm(n)
  integer ls(n)
  integer marker(n)
  integer nadj
  integer nblks
  real nonz(maxnon)
  integer nzsub(maxnon)
  integer perm(n)
  integer rchset(n)
  real rhs(n)
  integer xadj(n+1)
  integer xblk(maxblk+1)
  integer xenv(n+1)
  integer xls(n+1)
  integer xnonz(n+1)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  Use the 1WD method.'
  write ( *, * ) ' '
!
!  Initialize the permutation vectors.
!
  do i = 1, n
    perm(i) = i
    invprm(i) = i
  end do
!
!  Store the adjacency information.
!
  call adj_set_1 ( iadj, maxadj, nadj, n, xadj )

  write ( *, * ) ' '
  write ( *, * ) 'There are NADJ=',nadj,' adjacency entries.'
!
!  Display adjacency information.
!
  call adj_print ( iadj, maxadj, n, xadj )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Generate the 1WD ordering.
!
  call gen1wd ( n, xadj, iadj, nblks, xblk, perm, xls, ls )

  write ( *, * ) ' '
  write ( *, * ) 'Number of blocks is ',nblks
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )

  write ( *, * ) ' '
  write ( *, * ) '  The envelope size is ', ienv
!
!  Compute the inverse ordering.
!
  call perm_inverse ( n, perm, invprm )
!
!  Print orderings
!
  write ( *, * ) ' '
  write ( *, * ) '    I    Perm(I)   InvPerm(I)'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), invprm(i)
  end do
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Determine the envelope.
!
  call fnbenv ( xadj, iadj, perm, invprm, nblks, xblk, xenv, &
    ienv, marker, rchset, n )
!
!  Set RHS, DIAG, ENV.
!
  call setsy1 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, &
    n, nonz, nzsub, rhs, xenv, xnonz )
!
!  Factor the system.
!
  call ts_factor ( nblks, xblk, father, diag, xenv, env, xnonz, &
    nonz, nzsub, first, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST01 - Fatal error!'
    write ( *, * ) '  The matrix is not positive definite.'
    return
  end if
!
!  Solve the system.
!
  call ts_solve ( nblks, xblk, diag, xenv, env, xnonz, nonz, &
    nzsub, rhs, n )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, '(5g14.6)' ) ( rhs(i), i = 1, n )

  return
end
subroutine setsy1 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, nonz, &
  nzsub, rhs, xenv, xnonz )
!
!*******************************************************************************
!
!! SETSY1 stores the numerical values defining problem 1.
!
!
!  There is only one nonzero right hand side entry.
!  The matrix diagonal entries are all 2.
!  The nonzero offdiagonal entries are all -1.
!
  integer maxadj
  integer maxenv
  integer n
!
  real diag(n)
  real env(maxenv)
  integer i
  integer iadj(maxadj)
  integer invprm(n)
  integer isub
  integer j
  integer xadj(n+1)
  integer jsub
  real nonz(*)
  integer nzsub(*)
  real rhs(n)
  real value
  integer xenv(n+1)
  integer xnonz(n+1)
!
!  First zero out storage.
!
  rhs(1:n) = 0.0
  diag(1:n) = 0.0
  env(1:maxenv) = 0.0
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0

  call addrhs ( invprm, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1
 
      jsub = iadj(j)
      value = - 1.0

      call addrqt ( isub, jsub, value, invprm, diag, xenv, env, &
        xnonz, nonz, nzsub, n )

    end do

  end do

  return
end
subroutine adj_set_1 ( iadj, maxadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_SET_1 sets up the adjacency structure for problem 1.
!
  integer, parameter :: nonz = 19

  integer maxadj
  integer n
!
  integer iadj(maxadj)
  integer i
  integer, parameter, dimension ( nonz ) :: ilist = (/ &
    0,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/)
  integer, parameter, dimension ( nonz) :: jlist = (/ &
    0,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/)
  integer nadj
  integer xadj(n+1)
!
  do i = 1, nonz
    call adj_set ( iadj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests the ND method.
!
!
!  This example involves the following equations:
!
!  2*x1  - x10       = 0
!  2*x2  - x9  - x10 = 0
!  2*x3  - x8  - x9  = 0
!  2*x4  - x7  - x8  = 0
!  2*x5  - x6  - x7  = 0
!  2*x6  - x5        = 11
!  2*x7  - x4  - x5  = 0
!  2*x8  - x3  - x4  = 0
!  2*x9  - x2  - x3  = 0
!  2*x10 - x1  - x2  = 0
!
!  with solution (1,3,5,7,9,10,8,6,4,2).
!
  integer, parameter :: maxadj=300
  integer, parameter :: maxenv = 300
  integer, parameter :: maxnzsub = 300
  integer, parameter :: n = 10
!
  integer iadj(maxadj)
  real diag(n)
  real env(maxenv)
  integer first(n)
  integer i
  integer iband
  integer ienv
  integer invprm(n)
  integer link(n)
  integer lnz(99)
  integer ls(n)
  integer maxsub
  integer mrglnk(n)
  integer nadj
  integer nofnz
  integer nzsub(maxnzsub)
  integer perm(n)
  integer rchlnk(n)
  real rhs(n)
  integer xadj(n+1)
  integer xlnz(n+1)
  integer xls(n+1)
  integer xnzsub(99)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  Use the ND method.'
  write ( *, * ) ' '
!
!  Initialize the permutation vectors.
!
  do i = 1, n
    perm(i) = i
    invprm(i) = i
  end do
!
!  Store the adjacency information.
!
  call adj_set_2 ( iadj, maxadj, nadj, n, xadj )

  write ( *, * ) ' '
  write ( *, * ) 'There are NADJ=',nadj,' adjacency entries.'
!
!  Display adjacency information.
!
  call adj_print ( iadj, maxadj, n, xadj )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Generate the 1WD ordering.
!
  call gennd ( n, xadj, iadj, perm, xls, ls )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )

  write ( *, * ) ' '
  write ( *, * ) '  The envelope size is ', ienv
!
!  Compute the inverse ordering.
!
  call perm_inverse ( n, perm, invprm )
!
!  Print orderings
!
  write ( *, * ) ' '
  write ( *, * ) '    I    Perm(I)   InvPerm(I)'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), invprm(i)
  end do
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Symbolic factorization.
!
  maxsub = maxnzsub

  call smb_factor ( n, xadj, iadj, perm, invprm, xlnz, nofnz, &
    xnzsub, nzsub, maxsub, rchlnk, mrglnk )
!
!  Set RHS, DIAG, ENV.
!
  call setsy2 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, &
    n, rhs, nzsub, xnzsub, lnz, xlnz )
!
!  Factor the system.
!
  call gs_factor ( n, xlnz, lnz, xnzsub, nzsub, diag, link, first )
!
!  Solve the system.
!
  call gs_solve ( n, xlnz, lnz, xnzsub, nzsub, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy2 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, rhs, &
  nzsub, xnzsub, lnz, xlnz )
!
!*******************************************************************************
!
!! SETSY2 stores the numerical values defining problem 2.
!
!
!  There is only one nonzero right hand side entry.
!  The matrix diagonal entries are all 2.
!  The nonzero offdiagonal entries are all -1.
!
  integer maxadj
  integer maxenv
  integer n
!
  real diag(n)
  real env(maxenv)
  integer i
  integer iadj(maxadj)
  integer invprm(n)
  integer isub
  integer j
  integer jsub
  integer lnz(*)
  integer nzsub(*)
  real rhs(n)
  real value
  integer xadj(n+1)
  integer xlnz(*)
  integer xnzsub(n)
!
!  First zero out storage.
!
  rhs(1:n) = 0.0
  diag(1:n) = 0.0
  env(1:maxenv) = 0.0
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0

  call addrhs ( invprm, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1
 
      jsub = iadj(j)
      value = - 1.0

      call addcom ( isub, jsub, value, invprm, diag, lnz, xlnz, nzsub, xnzsub, &
        n )

    end do

  end do

  return
end
subroutine adj_set_2 ( iadj, maxadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_SET_2 sets up the adjacency structure for problem 2.
!
  integer, parameter :: nonz = 19

  integer maxadj
  integer n
!
  integer iadj(maxadj)
  integer i
  integer, parameter, dimension ( nonz ) :: ilist = (/ &
    0,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/)
  integer, parameter, dimension ( nonz ) :: jlist = (/ &
    0,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/)
  integer nadj
  integer xadj(n+1)
!
  do i = 1, nonz
    call adj_set ( iadj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests the QMD method.
!
!
!  This example involves the following equations:
!
!  2*x1  - x10       = 0
!  2*x2  - x9  - x10 = 0
!  2*x3  - x8  - x9  = 0
!  2*x4  - x7  - x8  = 0
!  2*x5  - x6  - x7  = 0
!  2*x6  - x5        = 11
!  2*x7  - x4  - x5  = 0
!  2*x8  - x3  - x4  = 0
!  2*x9  - x2  - x3  = 0
!  2*x10 - x1  - x2  = 0
!
!  with solution (1,3,5,7,9,10,8,6,4,2).
!
  integer, parameter :: maxadj = 300
  integer, parameter :: maxenv = 300
  integer, parameter :: maxnzsub = 300
  integer, parameter :: n = 10
!
  integer iadj(maxadj)
  integer iadj2(maxadj)
  real diag(n)
  real env(maxenv)
  integer first(n)
  integer i
  integer iband
  integer ienv
  integer invprm(n)
  integer link(n)
  integer lnz(99)
  integer marker(n)
  integer maxsub
  integer mrglnk(n)
  integer nadj
  integer nbrhd(n)
  integer nofnz
  integer nofsub
  integer nzsub(maxnzsub)
  integer perm(n)
  integer qlink(n)
  integer qsize(n)
  integer rchlnk(n)
  integer rchset(n)
  real rhs(n)
  integer xadj(n+1)
  integer xlnz(n+1)
  integer xnzsub(99)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  Use the QMD method.'
  write ( *, * ) ' '
!
!  Initialize the permutation vectors.
!
  do i = 1, n
    perm(i) = i
    invprm(i) = i
  end do
!
!  Store the adjacency information.
!
  call adj_set_3 ( iadj, maxadj, nadj, n, xadj )

  write ( *, * ) ' '
  write ( *, * ) 'There are NADJ=',nadj,' adjacency entries.'
!
!  Display adjacency information.
!
  call adj_print ( iadj, maxadj, n, xadj )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Generate the QMD ordering.
!
  iadj2 = iadj

  call genqmd ( n, xadj, iadj2, perm, invprm, marker, rchset, nbrhd, qsize, &
    qlink, nofsub )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )

  write ( *, * ) ' '
  write ( *, * ) '  The envelope size is ', ienv
!
!  Compute the inverse ordering.
!
  call perm_inverse ( n, perm, invprm )
!
!  Print orderings
!
  write ( *, * ) ' '
  write ( *, * ) '    I    Perm(I)   InvPerm(I)'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), invprm(i)
  end do
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Symbolic factorization.
!
  maxsub = maxnzsub

  call smb_factor ( n, xadj, iadj, perm, invprm, xlnz, nofnz, xnzsub, nzsub, &
    maxsub, rchlnk, mrglnk )
!
!  Set RHS, DIAG, ENV.
!
  call setsy3 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, rhs, nzsub, &
    xnzsub, lnz, xlnz )
!
!  Factor the system.
!
  call gs_factor ( n, xlnz, lnz, xnzsub, nzsub, diag, link, first )
!
!  Solve the system.
!
  call gs_solve ( n, xlnz, lnz, xnzsub, nzsub, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy3 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, rhs, &
  nzsub, xnzsub, lnz, xlnz )
!
!*******************************************************************************
!
!! SETSY3 stores the numerical values defining problem 3.
!
!
!  There is only one nonzero right hand side entry.
!  The matrix diagonal entries are all 2.
!  The nonzero offdiagonal entries are all -1.
!
  integer maxadj
  integer maxenv
  integer n
!
  real diag(n)
  real env(maxenv)
  integer i
  integer iadj(maxadj)
  integer invprm(n)
  integer isub
  integer j
  integer jsub
  integer lnz(*)
  integer nzsub(*)
  real rhs(n)
  real value
  integer xadj(n+1)
  integer xlnz(*)
  integer xnzsub(n)
!
!  First zero out storage.
!
  rhs(1:n) = 0.0
  diag(1:n) = 0.0
  env(1:maxenv) = 0.0
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0

  call addrhs ( invprm, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1
 
      jsub = iadj(j)
      value = - 1.0

      call addcom ( isub, jsub, value, invprm, diag, lnz, xlnz, nzsub, &
        xnzsub, n )

    end do

  end do

  return
end
subroutine adj_set_3 ( iadj, maxadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_SET_3 sets up the adjacency structure for problem 3.
!
  integer, parameter :: nonz = 19

  integer maxadj
  integer n
!
  integer iadj(maxadj)
  integer i
  integer ilist(nonz)
  integer jlist(nonz)
  integer nadj
  integer xadj(n+1)
!
  data ilist /0,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/
  data jlist /0,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/
!
  do i = 1, nonz
    call adj_set ( iadj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests the RCM method.
!
!
!  The first example involves the following equations:
!
!  2*x1  - x10       = 0
!  2*x2  - x9  - x10 = 0
!  2*x3  - x8  - x9  = 0
!  2*x4  - x7  - x8  = 0
!  2*x5  - x6  - x7  = 0
!  2*x6  - x5        = 11
!  2*x7  - x4  - x5  = 0
!  2*x8  - x3  - x4  = 0
!  2*x9  - x2  - x3  = 0
!  2*x10 - x1  - x2  = 0
!
!  with solution (1,3,5,7,9,10,8,6,4,2).
!
  integer, parameter :: maxadj = 3000
  integer, parameter :: maxenv = 3000
  integer, parameter :: n = 10
!
  integer bandw
  real diag(n)
  real env(maxenv)
  integer i
  integer iadj(maxadj)
  integer iband
  integer ienv
  integer ierror
  integer invprm(n)
  integer nadj
  integer perm(n)
  real rhs(n)
  integer xadj(n+1)
  integer xenv(n+1)
  integer xls(n+1)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  Use the RCM method.'
  write ( *, * ) ' '
!
!  Initialize the permutation vectors.
!
  do i = 1, n
    perm(i) = i
    invprm(i) = i
  end do
!
!  Store the adjacency information.
!
  call adj_set_4 ( iadj, maxadj, nadj, n, xadj )

  write ( *, * ) ' '
  write ( *, * ) 'There are NADJ=', nadj, ' adjacency entries.'
!
!  Display adjacency information
!
  call adj_print ( iadj, maxadj, n, xadj )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Generate the RCM ordering.
!
  call genrcm ( n, xadj, iadj, perm, xls )

  write ( *, * ) ' '
  write ( *, * ) 'The node used as starting point was originally'
  write ( *, * ) 'labeled ',perm(n)
!
!  Get inverse ordering
!
  call perm_inverse ( n, perm, invprm )
!
!  Print orderings
!
  write ( *, * ) ' '
  write ( *, * ) '    I    Perm(I)   InvPerm(I)'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), invprm(i)
  end do
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Compute the envelope.
!
  call fnenv ( n, xadj, iadj, perm, invprm, xenv, ienv, bandw )

  write ( *, * ) ' '
  write ( *, * ) '  The envelope size is ', ienv
  write ( *, * ) '  The bandwidth is ', bandw
!
!  Set RHS, DIAG, ENV.
!
  call setsy4 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, rhs, xenv )
!
!  Factor the matrix.
!
  call es_factor ( n, xenv, env, diag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST04 - Fatal error!'
    write ( *, * ) '  The matrix is not positive definite.'
    return
  end if
!
!  Solve the system.
!
  call el_solve ( n, xenv, env, diag, rhs )

  call eu_solve ( n, xenv, env, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy4 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, rhs, &
  xenv )
!
!*******************************************************************************
!
!! SETSY4 stores the numerical values defining problem 4.
!
!
!  There is only one nonzero right hand side entry.
!  The matrix diagonal entries are all 2.
!  The nonzero offdiagonal entries are all -1.
!
  integer maxadj
  integer maxenv
  integer n
!
  real diag(n)
  real env(maxenv)
  integer i
  integer iadj(maxadj)
  integer invprm(n)
  integer isub
  integer j
  integer xadj(n+1)
  integer jsub
  real rhs(n)
  real value
  integer xenv(n+1)
!
!  First zero out storage.
!
  rhs(1:n) = 0.0
  diag(1:n) = 0.0
  env(1:maxenv) = 0.0
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0

  call addrhs ( invprm, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = iadj(j)
      value = - 1.0

      call addrcm ( isub, jsub, value, invprm, diag, xenv, env, n )

    end do

  end do

  return
end
subroutine adj_set_4 ( iadj, maxadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_SET_4 sets up the adjacency structure for problem 4.
!
  integer, parameter :: nonz = 19

  integer maxadj
  integer n
!
  integer iadj(maxadj)
  integer i
  integer ilist(nonz)
  integer jlist(nonz)
  integer nadj
  integer xadj(n+1)
!
  data ilist /0,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/
  data jlist /0,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/
!
  do i = 1, nonz
    call adj_set ( iadj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests the RQT method.
!
!
!  This example involves the following equations:
!
!  2*x1  - x10       = 0
!  2*x2  - x9  - x10 = 0
!  2*x3  - x8  - x9  = 0
!  2*x4  - x7  - x8  = 0
!  2*x5  - x6  - x7  = 0
!  2*x6  - x5        = 11
!  2*x7  - x4  - x5  = 0
!  2*x8  - x3  - x4  = 0
!  2*x9  - x2  - x3  = 0
!  2*x10 - x1  - x2  = 0
!
!  with solution (1,3,5,7,9,10,8,6,4,2).
!
  integer, parameter :: maxadj = 300
  integer, parameter :: maxblk = 10
  integer, parameter :: maxenv = 300
  integer, parameter :: maxnon = 300
  integer, parameter :: n = 10
!
  integer iadj(maxadj)
  real diag(n)
  real env(maxenv)
  integer father(n)
  integer first(n)
  integer i
  integer iband
  integer ienv
  integer ierror
  integer invprm(n)
  integer ls(n)
  integer nadj
  integer nblks
  integer nodlvl(n)
  integer nofnz
  real nonz(maxnon)
  integer nzsub(maxnon)
  integer perm(n)
  real rhs(n)
  integer xadj(n+1)
  integer xblk(maxblk+1)
  integer xenv(n+1)
  integer xls(n+1)
  integer xnonz(n+1)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  Use the RQT method.'
  write ( *, * ) ' '
!
!  Initialize the permutation vectors.
!
  do i = 1, n
    perm(i) = i
    invprm(i) = i
  end do
!
!  Store the adjacency information.
!
  call adj_set_5 ( iadj, maxadj, nadj, n, xadj )
  write ( *, * ) ' '
  write ( *, * ) 'There are NADJ=', nadj, ' adjacency entries.'
!
!  Display adjacency information.
!
  call adj_print ( iadj, maxadj, n, xadj )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Generate the RQT ordering.
!
  call genrqt ( n, xadj, iadj, nblks, xblk, perm, xls, ls, nodlvl )

  write ( *, * ) ' '
  write ( *, * ) 'Number of blocks is ', nblks
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )

  write ( *, * ) ' '
  write ( *, * ) '  The envelope size is ', ienv

  call bshufl ( xadj, iadj, perm, nblks, xblk, xls, n )
!
!  Compute the inverse ordering.
!
  call perm_inverse ( n, perm, invprm )
!
!  Print orderings
!
  write ( *, * ) ' '
  write ( *, * ) '    I    Perm(I)   InvPerm(I)'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), invprm(i)
  end do
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Determine the quotient tree adjacency structure.
!
  call fntadj ( xadj, iadj, perm, nblks, xblk, father, n )
!
!  Determine the envelope index vector.
!
  call fntenv ( xadj, iadj, perm, invprm, nblks, xblk, xenv, ienv, n )

  nofnz = maxnon
  call fnofnz ( xadj, iadj, perm, invprm, nblks, xblk, xnonz, nzsub, nofnz, n )
!
!  Set RHS, DIAG, ENV.
!
  call setsy5 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, nonz, &
    nzsub, rhs, xenv, xnonz )
!
!  Factor the system.
!
  call ts_factor ( nblks, xblk, father, diag, xenv, env, xnonz, nonz, nzsub, &
    first, n, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST05 - Fatal error!'
    write ( *, * ) '  The matrix is not positive definite.'
    return
  end if
!
!  Solve the system.
!
  call ts_solve ( nblks, xblk, diag, xenv, env, xnonz, nonz, nzsub, rhs, n )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy5 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, nonz, &
  nzsub, rhs, xenv, xnonz )
!
!*******************************************************************************
!
!! SETSY5 stores the numerical values defining problem 5.
!
!
!  There is only one nonzero right hand side entry.
!  The matrix diagonal entries are all 2.
!  The nonzero offdiagonal entries are all -1.
!
  integer maxadj
  integer maxenv
  integer n
!
  real diag(n)
  real env(maxenv)
  integer i
  integer iadj(maxadj)
  integer invprm(n)
  integer isub
  integer j
  integer xadj(n+1)
  integer jsub
  real nonz(*)
  integer nzsub(*)
  real rhs(n)
  real value
  integer xenv(n+1)
  integer xnonz(n+1)
!
!  First zero out storage.
!
  rhs(1:n) = 0.0
  diag(1:n) = 0.0
  env(1:maxenv) = 0.0
!
!  Set the nonzero elements of the right hand side vector.
!
  isub = 6
  value = 11.0

  call addrhs ( invprm, isub, n, rhs, value )
!
!  Set the diagonal entries of the matrix.
!
  diag(1:n) = 2.0
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1
 
      jsub = iadj(j)
      value = - 1.0

      call addrqt ( isub, jsub, value, invprm, diag, xenv, env, xnonz, nonz, &
        nzsub, n )

    end do

  end do

  return
end
subroutine adj_set_5 ( iadj, maxadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_SET_3 sets up the adjacency structure for problem 5.
!
  integer, parameter :: nonz = 19

  integer maxadj
  integer n
!
  integer iadj(maxadj)
  integer i
  integer ilist(nonz)
  integer jlist(nonz)
  integer nadj
  integer xadj(n+1)
!
  data ilist /0,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/
  data jlist /0,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/
!
  do i = 1, nonz
    call adj_set ( iadj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 tests the RCM method.
!
!
!  This example corresponds to an implementation of the
!  5 point finite difference operator that approximates
!  the second derivative at a node by a centered difference.
!
!  The 8 by 8 grid is poorly numbered, to demonstrate how
!  SPARSPAK is able to improve a terrible bandwidth.
!
!  Here is the grid as defined by the user
!
!    18 02 57 41 11 29 48 21
!    49 38 12 25 45 63 37 01
!    30 62 53 07 60 52 15 56
!    08 46 64 16 34 28 24 40
!    42 22 26 33 17 06 44 10
!    58 13 05 61 27 51 59 32
!    03 35 54 43 23 14 36 47
!    19 50 31 09 39 55 04 20
!
!  The matrix consists of 4's on the diagonal, -1's on the
!  off-diagonals corresponding to the 4 immediate neighbors.
!  the right hand side is 2 on the corners, 1 on the sides,
!  and 0 in the interior.  The solution is (1,1,1,...,1,1,1).
!
  integer, parameter :: maxadj = 3000
  integer, parameter :: maxenv = 3000
  integer, parameter :: n = 64
!
  integer iadj(maxadj)
  integer bandw
  real diag(n)
  real env(maxenv)
  integer i
  integer iband
  integer ienv
  integer ierror
  integer invprm(n)
  integer nadj
  integer perm(n)
  real rhs(n)
  integer xadj(n+1)
  integer xenv(n+1)
  integer xls(n+1)
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  Use the RCM method.'
  write ( *, * ) ' '
!
!  Initialize the permutation vectors.
!
  do i = 1, n
    perm(i) = i
    invprm(i) = i
  end do
!
!  Store the adjacency information.
!
  call adj_set_6 ( iadj, maxadj, nadj, n, xadj )

  write ( *, * ) ' '
  write ( *, * ) 'There are NADJ=',nadj,' adjacency entries.'
!
!  Display the adjacency information
!
  call adj_print ( iadj, maxadj, n, xadj )
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Generate the RCM ordering.
!
  call genrcm ( n, xadj, iadj, perm, xls )

  write ( *, * ) ' '
  write ( *, * ) 'The node used as starting point was originally'
  write ( *, * ) 'labeled ', perm(n)
!
!  Get the inverse ordering.
!
  call perm_inverse ( n, perm, invprm )
!
!  Print the orderings.
!
  write ( *, * ) ' '
  write ( *, * ) '    I    Perm(I)   InvPrm(I)'
  write ( *, * ) ' '

  do i = 1, n
    write ( *, '(3i6)' ) i, perm(i), invprm(i)
  end do
!
!  Get a picture of the matrix.
!
  call shomat ( iadj, iband, ienv, invprm, maxadj, n, perm, xadj )
!
!  Compute the envelope.
!
  call fnenv ( n, xadj, iadj, perm, invprm, xenv, ienv, bandw )

  write ( *, * ) ' '
  write ( *, * ) '  The envelope size is ',ienv
  write ( *, * ) '  The bandwidth is ',bandw
!
!  Set the right hand side and the matrix.
!
  call setsy6 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, rhs, xenv )
!
!  Factor the matrix.
!
  call es_factor ( n, xenv, env, diag, ierror )

  if ( ierror /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'TEST02 - Fatal error!'
    write ( *, * ) '  The matrix is not positive definite.'
    return
  end if
!
!  Solve the system.
!
  call el_solve ( n, xenv, env, diag, rhs )

  call eu_solve ( n, xenv, env, diag, rhs )
!
!  Unpermute the solution.
!
  call perm_rv ( n, rhs, perm )

  write ( *, * ) ' '
  write ( *, * ) 'Solution:'
  write ( *, '(5g14.6)' ) rhs(1:n)

  return
end
subroutine setsy6 ( diag, env, iadj, invprm, xadj, maxadj, maxenv, n, rhs, &
  xenv )
!
!*******************************************************************************
!
!! SETSY6 stores the numerical values defining problem 6.
!
  integer maxadj
  integer maxenv
  integer n
!
  real diag(n)
  real env(maxenv)
  integer i
  integer iadj(maxadj)
  integer, parameter, dimension ( 28 ) :: ilist = (/ &
    1, 2, 3, 4, 8, 9,10,11,18,19,20,21,29,30,31,32,39, &
   40,41,42,47,48,49,50,55,56,57,58 /)
  integer invprm(n)
  integer isub
  integer j
  integer xadj(n+1)
  integer jsub
  real rhs(n)
  real, parameter, dimension ( 28 ) :: rlist = (/ &
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,1.0,1.0,1.0,1.0,1.0, &
    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0 /)
  real value
  integer xenv(n+1)
!
!  Zero out the arrays.
!
  rhs(1:n) = 0.0
  diag(1:n) = 4.0
  env(1:maxenv) = 0.0
!
!  Set the nonzero elements of the right hand side.
!
  do i = 1, 28
    call addrhs ( invprm, ilist(i), n, rhs, rlist(i) )
  end do
!
!  Set the off diagonal terms of the matrix.
!
  do i = 1, n

    isub = i

    do j = xadj(i), xadj(i+1) - 1

      jsub = iadj(j)
      value = - 1.0

      call addrcm ( isub, jsub, value, invprm, diag, xenv, env, n )

    end do

  end do

  return
end
subroutine adj_set_6 ( iadj, maxadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_SET_6 sets up the adjacency structure for problem 6.
!
  integer maxadj
  integer n
!
  integer iadj(maxadj)
  integer i
  integer, parameter, dimension ( 113 ) :: ilist = (/ &
     0,13,16,17,18,19,20,21,22,23,24,25,25,26,26,27,27,28,28,29, &
    30,31,32,33,33,33,34,34,34,35,35,36,36,37,37,38,38,39,39,40, &
    40,41,41,42,42,43,43,44,44,44,45,45,46,46,47,47,47,48,48,48, &
    49,49,49,50,50,50,51,51,51,52,52,53,53,54,54,54,54,55,55,55, &
    56,56,56,57,57,57,58,58,58,59,59,59,59,60,60,60,60,61,61,61, &
    61,62,62,62,62,63,63,63,63,64,64,64,64 /)
  integer xadj(n+1)
  integer, parameter, dimension ( 113 ) :: jlist = (/ &
     0, 5, 7, 6, 2, 3, 4, 1,13,14,15, 7,12, 5,22,17,23, 6,24,11, &
     8, 9,10,16,17,26,16,17,28, 3,13, 4,14, 1,15, 2,12, 9,23,10, &
    24,11,25, 8,22, 9,23, 6,10,24,11,25, 8,22,20,32,36,21,29,37, &
    18,30,38,19,31,35, 6,14,27,15,28, 7,12, 5,31,35,43, 4,14,39, &
     1,15,40, 2,12,41, 3,13,42,32,36,44,51, 7,34,45,52, 5,27,33, &
    43,30,38,46,53,29,37,45,52,16,26,46,53 /)
  integer nadj
!
  do i = 1, 113
    call adj_set ( iadj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
  end do

  return
end
