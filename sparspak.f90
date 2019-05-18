!  sparspak.f90  17 June 2000
!
subroutine addcom ( isub, jsub, value, invprm, diag, xlnz, ixlnz, nzsub, &
  xnzsub, n )
!
!*******************************************************************************
!
!! ADDCOM can add values to a matrix stored in compressed storage scheme.
!
!
!  Discussion:
!
!    The routine is called once the equations and variables have
!    been reordered, to add numbers to the matrix.  The storage area 
!    should be zeroed out before this routine is called.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ISUB, JSUB, the row and column of the matrix to
!    which the value is to be added.  Note that the matrix is assumed
!    to be symmetric, and only the lower triangle is stored.  Hence,
!    if ISUB < JSUB, the routine ignores the request.
!
!    Input, real VALUE, the quantity to be added to A(ISUB,JSUB).
!
!    Input, integer INVPRM(N), the inverse ordering, which should have
!    been created by calling INVRSE once the reordering is set.
!
!    Input/output, real DIAG(N), the diagonal elements of the matrix.
!
!    Input/output, real XLNZ(*), the nonzero subdiagonal elements of the matrix.
!
!    Input, integer IXLNZ(N+1), NZSUB(*), XNZSUB(N), data structures which
!    define the compressed storage scheme.
!
!    Input, integer N, the order of the matrix.
!
  integer n
!
  real diag(n)
  integer i
  integer invprm(n)
  integer isub
  integer ixlnz(n+1)
  integer j
  integer jsub
  integer k
  integer kstop
  integer kstrt
  integer ksub
  integer nzsub(*)
  real value
  real xlnz(*)
  integer xnzsub(n)
!
!  Figure out the current locations of the given row and column.
!
  i = invprm(isub)
  j = invprm(jsub)
!
!  If the entry is on the diagonal, update DIAG.
!
  if ( i == j ) then
    diag(i) = diag(i) + value
    return
  end if
!
!  If the entry is above the diagonal, don't store it at all.
!
  if ( i < j ) then
    return
  end if

  kstrt = ixlnz(j)
  kstop = ixlnz(j+1)-1

  if ( kstop < kstrt ) then
    write ( *, * ) ' '
    write ( *, * ) 'ADDCOM - Fatal error!'
    write ( *, * ) '  The IXLNZ array is incorrect.'
    write ( *, * ) '  IXLNZ(J) =     ', ixlnz(j)
    write ( *, * ) '  IXLNZ(J+1)-1 = ', ixlnz(j+1)-1
    stop
  end if

  ksub = xnzsub(j)

  do k = kstrt, kstop

    if ( nzsub(ksub) == i ) then
      xlnz(k) = xlnz(k) + value
      return
    end if

    ksub = ksub + 1

  end do

  write ( *, * ) ' '
  write ( *, * ) 'ADDCOM - Fatal error!'
  write ( *, * ) '  No storage was set aside for entry ISUB, JSUB'
  write ( *, * ) '  where ISUB = ', isub, ' and JSUB = ', jsub

  stop
end
subroutine addrcm ( isub, jsub, value, invprm, diag, xenv, env, n )
!
!*******************************************************************************
!
!! ADDRCM can add values to a matrix stored in the RCM scheme.
!
!
!  Discussion:
!
!    Since this routine only adds VALUE to the current matrix entry, the 
!    matrix should be zeroed out first.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ISUB, JSUB, the row and column of the matrix to
!    which the value is to be added.  Note that the matrix is assumed
!    to be symmetric, and only the lower triangle is stored.  Hence,
!    if ISUB < JSUB, the routine ignores the request.
!
!    Input, real VALUE, the number to be added.
!
!    Input, integer INVPRM(N), the inverse variable ordering.
!
!    Input/output, real DIAG(N), the diagonal of the matrix.
!
!    Input, integer XENV(N+1), describes the envelope structure.
!
!    Input/output, real ENV(*), contains the nonzeros of the matrix.
!
!    Input, integer N, the order of the matrix.
!
  integer n
!
  real diag(n)
  real env(*)
  integer i
  integer invprm(n)
  integer isub
  integer j
  integer jsub
  integer k
  real value
  integer xenv(n+1)
!
  i = invprm(isub)
  j = invprm(jsub)

  if ( i < j ) then
    return
  end if

  if ( i == j ) then
    diag(i) = diag(i) + value
    return
  end if

  k = xenv(i+1) - i + j

  if ( k < xenv(i) ) then
    write ( *, * ) ' '
    write ( *, * ) 'ADDRCM - Fatal error!'
    write ( *, * ) '  Indices outside of envelope.'
    write ( *, * ) '  ISUB = ',isub,' JSUB = ',jsub
    stop
  end if

  env(k) = env(k) + value

  return
end
subroutine addrhs ( invprm, isub, n, rhs, value )
!
!*******************************************************************************
!
!! ADDRHS adds a quantity to a specific entry of the right hand side.
!
!
!  Discussion:
!
!    After the equations have been reordered, it may be desired to
!    alter one of the entries of the right hand side.  This routine
!    carries out this operation after adjusting for the reordering.
!
!  Modified:
!
!    29 April 2000
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer INVPRM(N), the inverse permutation, constructed
!    by PERM_INVERSE.
!
!    Input, integer ISUB, the index of the entry of RHS to be modified.
!
!    Input, integer N, the order of the system.
!
!    Input/output, real RHS(N), the right hand side.  This vector
!    has been permuted.  
!
!    Input, real VALUE, the quantity to be added to the right hand side.
!
  integer n
!
  integer i
  integer invprm(n)
  integer isub
  real rhs(n)
  real value
!
  i = invprm(isub)

  if ( 1 <= i .and. i <= n ) then
    rhs(i) = rhs(i) + value
  end if

  return
end
subroutine addrqt ( isub, jsub, value, invprm, diag, xenv, env, xnonz, nonz, &
  nzsubs, n )
!
!*******************************************************************************
!
!! ADDRQT can add values to a matrix stored in the implicit block storage scheme.
!
!
!  Discussion:
!
!    Since the routine only adds new values to those currently in storage, the
!    space used to store the matrix must be initialized to 0 before numerical 
!    values are supplied.  This routine can be used with the RQT and 1WD 
!    methods.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ISUB, JSUB, the row and column of the matrix to
!    which the value is to be added.  Note that the matrix is assumed
!    to be symmetric, and only the lower triangle is stored.  Hence,
!    if ISUB < JSUB, the routine ignores the request.
!
!    Input, real VALUE, the number to be added.
!
!    Input, integer INVPRM(N), the inverse variable ordering.
!
!    Input/output, real DIAG(N), the diagonal of the matrix.
!
!    Input, integer XENV(N+1), describes the envelope structure of the diagonal
!    blocks.
!
!    Input/output, real ENV(*), contains the nonzeros of the matrix.
!
!    Input/output, integer XNONZ(N+1), real NONZ(*), integer NZSUBS(*),
!    levels structure containing the off-block diagonal parts of the rows
!    of the lower triangle of the original matrix.
!
!    Input, integer N, the order of the matrix.
!
  integer n
!
  real diag(n)
  real env(*)
  integer i
  integer invprm(n)
  integer isub
  integer j
  integer jsub
  integer k
  integer kstop
  integer kstrt
  real nonz(*)
  integer nzsubs(*)
  real value
  integer xenv(n+1)
  integer xnonz(n+1)
!
  i = invprm(isub)
  j = invprm(jsub)
!
!  Ignore superdiagonal entries.
!
  if ( i < j ) then
    return
  end if
!
!  Diagonal entries.
!
  if ( i == j ) then
    diag(i) = diag(i) + value
    return
  end if
!
!  Entries within the diagonal envelope.
!
  k = xenv(i+1) - i + j

  if ( k >= xenv(i) ) then
    env(k) = env(k) + value
    return
  end if
!
!  The value goes outside the diagonal blocks.
!
  kstrt = xnonz(i)
  kstop = xnonz(i+1) - 1

  do k = kstrt, kstop

    if ( nzsubs(k) == j ) then
      nonz(k) = nonz(k) + value
      return
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'ADDRQT - Fatal error!'
  write ( *, * ) '  Lack of storage!'
  write ( *, * ) '  ISUB = ', isub
  write ( *, * ) '  JSUB = ', jsub

  stop
end
subroutine adj_print ( iadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_PRINT prints the adjacency information stored in XADJ, ADJNCY.  
!
!
!  Discussion:
!
!    The list has the form:
!
!    Row   Nonzeros
!
!    1       2   5   9
!    2       7   8   9   15   78   79   81  86  91  99
!          100 103
!    3      48  49  53
!
!
!  Parameters:
!
!    Input, integer IADJ(NADJ), the adjacency structure, which contains, 
!    for each row, the column indices of the nonzero entries.
!
!    Input, integer NADJ, the dimension of ADJNCY.
!
!    Input, integer N, the number of equations.
!
!    Input, integer XADJ(N+1), organizes the entries of ADJNCY
!    into rows.  The entries for row I are in entries XADJ(I)
!    through XADJ(I+1)-1.
!
  integer nadj
  integer n
!
  integer iadj(nadj)
  integer i
  integer j
  integer jhi
  integer jlo
  integer jmax
  integer jmin
  integer xadj(n+1)
!
  write ( *, * ) ' '
  write ( *, * ) 'ADJ_PRINT'
  write ( *, * ) '  Show adjacency structure of sparse matrix.'
  write ( *, * ) '  There are a total of ',nadj,' entries.'
  write ( *, * ) ' '
  write ( *, * ) 'Row         Nonzeros '
  write ( *, * ) ' '

  do i = 1, n
    jmin = xadj(i)
    jmax = xadj(i+1)-1
    do jlo = jmin, jmax, 10
      jhi = min ( jlo+9, jmax )
      if ( jlo == jmin ) then
        write ( *, '(i6,6x,10i6)' ) i, ( iadj(j), j = jlo, jhi )
      else
        write ( *, '(6x,6x,10i6)' )  ( iadj(j), j = jlo, jhi )
      end if
    end do
  end do

  return
end
subroutine adj_set ( iadj, irow, jcol, maxadj, nadj, n, xadj )
!
!*******************************************************************************
!
!! ADJ_SET sets up the adjacency information XADJ and ADJNCY.
!
!
!  Discussion:
!
!    The routine records the locations of each nonzero element,
!    one at a time.
!
!    The first call  for a given problem should be with
!    irow = jcol = 0, which is a flag to clear out
!    the old information, and set xadj(i) = 1 for i = 1,n+1.
!
!    Thereafter, call with values of irow and jcol
!    for each entry a(irow,jcol) of the matrix that is nonzero.
!
!    Diagonal entries (where I = J) are not stored in ADJNCY, and so
!    you do not need to report them.
!
!    The matrix is assumed to be symmetric, and so you only have to give
!    the subdiagonal or superdiagonal entries.
!
!    repeated calls with the same values of irow and jcol do not
!    hurt.  no extra storage will be allocated.
!
!  Modified:
!
!    29 April 2000
!
!  Parameters:
!
!    Input/output, integer IADJ(NADJ), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer IADJ(MAXADJ), stores adjacency information.
!
!    Input, integer IROW, JCOL, the row and column indices of a nonzero
!    entry of the matrix.  The special case IROW = JCOL = 0 signals
!    initialization.  If IROW = JCOL, a diagonal entry is being specified,
!    but this is unnecessary since space is always set aside for diagonal
!    entries.  The matrix is assumed to be symmetric, so you only need
!    to specify the subdiagonal nonzeros.
!
!    Input, integer MAXADJ, the dimension of IADJ.
!
!    Input/output, integer NADJ, the number of adjaceny entries in ADJNCY.
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
  integer maxadj
  integer n
!
  integer iadj(maxadj)
  integer i
  integer irow
  integer j
  integer jcol
  integer k
  integer kback
  integer khi
  integer klo
  integer nadj
  integer xadj(n+1)
!
  if ( irow == 0 .or. jcol == 0 ) then

    write ( *, * ) ' '
    write ( *, * ) 'ADJ_SET - Note:'
    write ( *, * ) '  Initializing adjacency information.'
    write ( *, * ) '  Number of equations N = ', n
    write ( *, * ) '  Maximum adjacency MAXADJ = ', maxadj

    nadj = 0
    xadj(1:n+1) = 1
    iadj(1:maxadj) = 0

    return

  end if

  if ( irow == jcol ) then
    return
  end if

  if ( xadj(n+1) > maxadj + 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ADJ_SET - Fatal error!'
    write ( *, * ) '  All available storage has been used.'
    write ( *, * ) '  No more information can be stored!'
    write ( *, * ) '  This error occurred for '
    write ( *, * ) '  Row IROW = ',irow
    write ( *, * ) '  Column JCOL = ',jcol
    stop
  end if

  if ( irow > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'ADJ_SET - Fatal error!'
    write ( *, * ) '  IROW > N.'
    write ( *, * ) '  IROW = ', irow
    write ( *, * ) '  N = ', n
    stop
  else if ( irow < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ADJ_SET - Fatal error!'
    write ( *, * ) '  IROW < 1.'
    write ( *, * ) '  IROW = ', irow
    stop
  else if ( jcol > n ) then
    write ( *, * ) ' '
    write ( *, * ) 'ADJ_SET - Fatal error!'
    write ( *, * ) '  JCOL > N.'
    write ( *, * ) '  JCOL = ', jcol
    write ( *, * ) '  N = ', n
    stop
  else if ( jcol < 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ADJ_SET - Fatal error!'
    write ( *, * ) '  JCOL < 1.'
    write ( *, * ) '  JCOL = ', jcol
    stop
  end if

  i = irow
  j = jcol

   20 continue
!
!  Search the adjacency entries already stored for row I,
!  to see if J has already been stored.
!
  klo = xadj(i)
  khi = xadj(i+1)-1

  do k = klo, khi

    if ( iadj(k) == j ) then

      if ( i == irow ) then
        i = jcol
        j = irow
        go to 20
      end if

      return

    end if

  end do
!
!  A new adjacency entry must be made.
!  Shift the later ADJNCY entries up one.
!
  do k = xadj(i+1), xadj(n+1)
    kback = xadj(n+1) + xadj(i+1) - k
    iadj(kback+1) = iadj(kback)
  end do
!
!  Insert the new entry.
!
  iadj(xadj(i+1)) = j
!
!  Update the XADJ pointers.
!
  do k = i+1, n+1
    xadj(k) = xadj(k) + 1
  end do

  nadj = xadj(n+1) - 1

  if ( i == irow ) then
    i = jcol
    j = irow
    go to 20
  end if

  return
end
subroutine bshufl ( xadj, iadj, perm, nblks, xblk, xls, n )
!
!*******************************************************************************
!
!! BSHUFL renumbers the nodes of each block to reduce its envelope.  
!
!
!  Discussion:
!
!    Nodes in a block with no neighbors in previous blocks are renumbered 
!    by SUBRCM before the others.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  (nblks,xblk)-the tree partitioning.
!
!  updated parameter
!  perm-the permutation vector. on return, it contains
!       the new permutation.
!
!  working vectors-
!
!  xls-index vector to a level structure.
!
  integer nblks
  integer n
!
  integer iadj(*)
  integer bnum(n)
  integer i
  integer ip
  integer istop
  integer istrt
  integer j
  integer jstop
  integer jstrt
  integer k
  integer mask(n)
  integer nabor
  integer nbrblk
  integer node
  integer nsubg
  integer perm(n)
  integer subg(n)
  integer xadj(n+1)
  integer xblk(nblks+1)
  integer xls(n+1)
!
  if ( nblks <= 0 ) then
    return
  end if
!
!  Find the block number for each variable and initialize MASK.
!
!  BNUM stores the block number of each variable.
!
!  MASK is used to prescribe a subgraph.
!
  do k = 1, nblks

    istrt = xblk(k)
    istop = xblk(k+1)-1

    do i = istrt, istop
      node = perm(i)
      bnum(node) = k
      mask(node) = 0
    end do

  end do
!
!  For each block, find those nodes with no neighbors
!  in previous blocks and accumulate them in SUBG.
!  They will be renumbered before others in the block.
!
  do k = 1, nblks

    istrt = xblk(k)
    istop = xblk(k+1)-1
    nsubg = 0

    do i = istrt, istop

      node = perm(i)
      jstrt = xadj(node)
      jstop = xadj(node+1)-1

      if ( jstop >= jstrt ) then

        do j = jstrt, jstop
          nabor = iadj(j)
          nbrblk = bnum(nabor)
          if ( nbrblk < k ) then
            go to 40
          end if
        end do

        nsubg = nsubg+1
        subg(nsubg) = node
        ip = istrt+nsubg-1
        perm(i) = perm(ip)

      end if

40        continue

    end do
!
!  SUBRCM renumbers the subgraph stored in (NSUBG, SUBG).
!
    if ( nsubg > 0 ) then
      call subrcm ( xadj, iadj, mask, nsubg, subg, perm(istrt), xls, n )
    end if

  end do

  return
end
subroutine degree ( root, xadj, iadj, mask, deg, iccsze, ls, n )
!
!*******************************************************************************
!
!! DEGREE computes the degrees of the nodes in the connected component.
!
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, is the node that defines the component.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  mask-specifies a section subgraph.
!
!  output parameters-
!  deg-array containing the degrees of the nodes in
!  the component.
!  iccsze-size of the component specifed by mask and root
!
!  working parameter-
!  ls-a temporary vector used to store the nodes of the
!  component level by level.
!
  integer n
!
  integer iadj(*)
  integer deg(n)
  integer i
  integer iccsze
  integer ideg
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer ls(n)
  integer lvlend
  integer lvsize
  integer mask(n)
  integer nbr
  integer node
  integer root
  integer xadj(n+1)
!
!  The array XADJ is used as a temporary marker to
!  indicate which nodes have been considered so far.
!
  ls(1) = root
  xadj(root) = -xadj(root)
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
   10 continue

  lbegin = lvlend+1
  lvlend = iccsze
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
!
  do i = lbegin, lvlend

    node = ls(i)
    jstrt = -xadj(node)
    jstop = abs ( xadj(node+1) ) - 1
    ideg = 0

    do j = jstrt, jstop

      nbr = iadj(j)

      if ( mask(nbr) /= 0 ) then

        ideg = ideg+1

        if ( xadj(nbr) >= 0 ) then
          xadj(nbr) = -xadj(nbr)
          iccsze = iccsze+1
          ls(iccsze) = nbr
        end if

      end if

    end do

    deg(node) = ideg

  end do
!
!  Compute the current level width.
!
  lvsize = iccsze - lvlend
!
!  If the current level width is nonzero, generate another level.
!
  if ( lvsize > 0 ) then
    go to 10
  end if
!
!  Reset XADJ to its correct sign and return.
!
  do i = 1, iccsze
    node = ls(i)
    xadj(node) = -xadj(node)
  end do

  return
end
subroutine el_solve ( n, xenv, env, diag, rhs )
!
!*******************************************************************************
!
!! EL_SOLVE solves a lower triangular system stored in the envelope format.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the number of equations.
!
!    Input, (xenv,env)-array pair for the envelope of L.
!
!    Input, real DIAG(N), the diagonal of the matrix.
!
!    Input/output, real RHS(N).  On input, the right hand side.
!    On output, the solution.
!
  integer n
!
  real diag(n)
  real env(*)
  integer i
  integer iband
  integer ifirst
  integer k
  integer kstop
  integer kstrt
  integer l
  integer last
  real rhs(n)
  real s
  integer xenv(n+1)
!
!  Find the position of the first nonzero in RHS and put it in IFIRST.
!
  ifirst = 0

   10 continue

  ifirst = ifirst + 1

  if ( rhs(ifirst) ==  0.0 ) then
    if ( ifirst >= n ) then
      return
    end if
    go to 10
  end if
!
!  LAST contains the position of the most recently
!  computed nonzero component of the solution.
!
  last = 0

  do i = ifirst, n

    iband = xenv(i+1)-xenv(i)
    iband = min ( iband, i-1 )

    s = rhs(i)
    l = i-iband
    rhs(i) = 0.0
!
!  row of the envelope is empty, or corresponding
!  components of the solution are all zeros.
!
    if ( iband /= 0 .and. last >= l ) then

      kstrt = xenv(i+1)-iband
      kstop = xenv(i+1)-1

      do k = kstrt, kstop
        s = s - env(k) * rhs(l)
        l = l + 1
      end do

    end if

    if ( s /= 0.0 ) then
      rhs(i) = s / diag(i)
      last = i
    end if

  end do

  return
end
subroutine es_factor ( n, xenv, env, diag, ierror )
!
!*******************************************************************************
!
!! ES_FACTOR factors a positive definite envelope matrix into l*l(transpose).
!
!
!  Discussion:
!
!    The algorithm used is the standard bordering method.
!
!  Modified:
!
!    29 April 2000
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    xenv-the envelope index vector.
!
!    updated parameters-
!    env-the envelope of l overwrites that of a.
!
!    diag-the diagonal of l overwrites that of a.
!
!    Output, integer IERROR, error flag.
!    0, no error, the factorization was carried out.
!    1, the matrix is not positive definite.
!
  integer n
!
  real diag(n)
  real env(*)
  integer i
  integer iband
  integer ierror
  integer ifirst
  integer ixenv
  integer j
  integer jstop
  real s
  real temp
  integer xenv(n+1)
!
  if ( diag(1) <= 0.0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ES_FACTOR - Fatal error!'
    write ( *, * ) '  The matrix is not positive definite.'
    ierror = 1
    return
  end if

  diag(1) = sqrt ( diag(1) )

  if ( n == 1 ) then
    ierror = 0
    return
  end if
!
!  Loop over rows 2,3,...,N of the matrix
!
  do i = 2, n

    ixenv = xenv(i)
    iband = xenv(i+1) - ixenv
!
!  Compute row I of the triangular factor.
!
    temp = diag(i)

    if ( iband /= 0 ) then

      ifirst = i - iband

      call el_solve ( iband, xenv(ifirst), env, diag(ifirst), env(ixenv) )

      jstop = xenv(i+1) - 1

      do j = ixenv, jstop
        s = env(j)
        temp = temp - s * s
      end do

    end if

    if ( temp <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'ES_FACTOR - Fatal error!'
      write ( *, * ) '  The matrix is not positive definite.'
      ierror = 1
      return
    end if

    diag(i) = sqrt ( temp )

  end do

  ierror = 0

  return
end
subroutine eu_solve ( n, xenv, env, diag, rhs )
!
!*******************************************************************************
!
!! EU_SOLVE solves an upper triangular system stored in the envelope format.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!  (xenv,env)-array pair for the envelope of u.
!  diag-array for the diagonal of u.
!
!  updated parameter-
!  rhs-on input, it contains the right hand side.
!  on output, it contains the solution vector.
!
  integer n
!
  real diag(n)
  real env(*)
  integer i
  integer iband
  integer k
  integer kstop
  integer kstrt
  integer l
  real rhs(n)
  real s
  integer xenv(n+1)
!
  i = n + 1

   10 continue

  i = i - 1

  if ( i == 0 ) then
    return
  end if

  if ( rhs(i) == 0.0 ) then
    go to 10
  end if

  s = rhs(i)/diag(i)
  rhs(i) = s
  iband = xenv(i+1)-xenv(i)

  if ( iband >= i ) then
    iband = i-1
  end if

  if ( iband == 0 ) then
    go to 10
  end if

  kstrt = i-iband
  kstop = i-1
  l = xenv(i+1)-iband

  do k = kstrt, kstop
    rhs(k) = rhs(k) - s * env(l)
    l = l+1
  end do

  go to 10
end
subroutine fnbenv ( xadj, iadj, perm, invprm, nblks, xblk, xenv, ienv, marker, &
  rchset, n )
!
!*******************************************************************************
!
!! FNBENV finds the exact envelope structure of the diagonal blocks of
!  the Cholesky factor of a permuted partitioned matrix.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  (perm,invprm)-the permutation vector and its inverse.
!  (nblks,xblk)-the partitioning.
!
!  output parameters -
!  xenv-the envelope index vector.
!  ienv-the size of the envelope.
!
!  working parameters-
!  mask-marks nodes that have been considered.
!  marker-is used by routine reach.
!  rchset-is used by the subroutine reach.
!  stores both reachable and neighborhood sets.
!
  integer nblks
  integer n
!
  integer iadj(*)
  integer blkbeg
  integer blkend
  integer i
  integer ienv
  integer ifirst
  integer inhd
  integer invprm(n)
  integer k
  integer marker(n)
  integer mask(n)
  integer newnhd
  integer nhdsze
  integer node
  integer perm(n)
  integer rchset(*)
  integer rchsze
  integer xadj(n+1)
  integer xblk(nblks+1)
  integer xenv(n+1)
!
  n = xblk(nblks+1)-1
  ienv = 1
!
!  MASK will mark nodes that have been considered.
!
  mask(1:n) = 0
  marker(1:n) = 1
!
!  Loop over the blocks.
!
  do k = 1, nblks

    nhdsze = 0
    blkbeg = xblk(k)
    blkend = xblk(k+1)-1

    do i = blkbeg, blkend
      node = perm(i)
      marker(node) = 0
    end do
!
!  Loop through the nodes in the current block.
!
    do i = blkbeg, blkend
      node = perm(i)
      call reach ( node, xadj, iadj, mask, marker, rchsze, rchset(blkbeg), &
        newnhd, rchset(nhdsze+1), n )
      nhdsze = nhdsze+newnhd
      ifirst = marker(node)
      ifirst = invprm(ifirst)
      xenv(i) = ienv
      ienv = ienv+i-ifirst
    end do
!
!  Reset marker values of nodes in nbrhd set.
!
    do inhd = 1, nhdsze
      node = rchset(inhd)
      marker(node) = 0
    end do
!
!  Reset marker and mask values of nodes in the current block.
!
    do i = blkbeg, blkend
      node = perm(i)
      marker(node) = 0
      mask(node) = 1
    end do

  end do

  xenv(n+1) = ienv
  ienv = ienv-1

  return
end
subroutine fndsep ( root, xadj, iadj, mask, nsep, sep, xls, ls, n )
!
!*******************************************************************************
!
!! FNDSEP finds a small separator for a connected component in a graph.
!
!
!  Discussion:
!
!    The connected component is specified by the MASK array.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that determines the masked component.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameters-
!  nsep-number of variables in the separator.
!  sep-vector containing the separator nodes.
!
!  updated parameter-
!  mask-nodes in the separator have their mask
!  values set to zero.
!
!  working parameters-
!  (xls,ls)-level structure pair for level structure
!  found by fnroot.
!
  integer n
!
  integer iadj(*)
  integer i
  integer j
  integer jstop
  integer jstrt
  integer ls(n)
  integer mask(n)
  integer midbeg
  integer midend
  integer midlvl
  integer mp1beg
  integer mp1end
  integer nbr
  integer nlvl
  integer node
  integer nsep
  integer root
  integer sep(*)
  integer xadj(n+1)
  integer xls(n+1)
!
  call fnroot ( root, xadj, iadj, mask, nlvl, xls, ls, n )
!
!  If the number of levels is less than 3, return the whole component
!  as the separator.
!
  if ( nlvl < 3 ) then

    nsep = xls(nlvl+1)-1

    do i = 1, nsep
      node = ls(i)
      sep(i) = node
      mask(node) = 0
    end do

    return

  end if
!
!  Find the middle level of the rooted level structure.
!
  midlvl = (nlvl+2)/2
  midbeg = xls(midlvl)
  mp1beg = xls(midlvl+1)
  midend = mp1beg-1
  mp1end = xls(midlvl+2)-1
!
!  The separator is obtained by including only those
!  middle-level nodes with neighbors in the middle+1
!  level. xadj is used temporarily to mark those
!  nodes in the middle+1 level.
!
  do i = mp1beg, mp1end
    node = ls(i)
    xadj(node) = -xadj(node)
  end do

  nsep = 0

  do i = midbeg, midend

    node = ls(i)
    jstrt = xadj(node)
    jstop = abs ( xadj(node+1) ) - 1

    do j = jstrt, jstop

      nbr = iadj(j)

      if ( xadj(nbr) <= 0 ) then
        nsep = nsep+1
        sep(nsep) = node
        mask(node) = 0
        go to 10
      end if

    end do

10      continue

  end do
!
!  Reset XADJ to its correct sign.
!
  do i = mp1beg, mp1end
    node = ls(i)
    xadj(node) = -xadj(node)
  end do

  return
end
subroutine fnenv ( n, xadj, iadj, perm, invprm, xenv, ienv, iband )
!
!*******************************************************************************
!
!! FNENV finds the envelope structure of a permuted matrix.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input/output, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  perm,invprm-arrays containing permutation data about
!  the reordered matrix.
!
!  output parameters-
!  xenv-index vector for the level structure
!  to be used to store the lower (or upper)
!  envelope of the reordered matrix.
!  ienv-is equal to xenv(n+1)-1.
!  iband-bandwidth of the reordered matrix.
!
  integer n
!
  integer iadj(*)
  integer i
  integer iband
  integer ienv
  integer ifirst
  integer invprm(n)
  integer iperm
  integer j
  integer jband
  integer jstop
  integer jstrt
  integer nabor
  integer perm(n)
  integer xadj(n+1)
  integer xenv(n+1)
!
  iband = 0
  ienv = 1

  do i = 1, n

    xenv(i) = ienv
    iperm = perm(i)
    jstrt = xadj(iperm)
    jstop = xadj(iperm+1)-1
!
!  Find the first nonzero in row I.
!
    if ( jstrt <= jstop ) then

      ifirst = i

      do j = jstrt, jstop

        nabor = iadj(j)
        nabor = invprm(nabor)

        if ( nabor < ifirst ) then
          ifirst = nabor
        end if

      end do

      jband = i - ifirst
      ienv = ienv + jband
      iband = max ( iband, jband )

    end if

  end do

  xenv(n+1) = ienv
  ienv = ienv-1

  return
end
subroutine fnlvls ( root, xadj, iadj, nodlvl, nlvl, xls, ls, n )
!
!*******************************************************************************
!
!! FNLVLS generates a rooted level structure for a masked connected
!  subgraph, rooted at a pseudo-peripheral node.  the level numbers
!  are recorded.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input/output, integer ROOT.  On input, with the array NODLVL, ROOT
!    specifies the component whose pseudo-peripheral node is to be found.
!    On output, ROOT contains the value of that pseudo-peripheral node.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameters-
!  nlvl-number of levels in the level structure found.
!  (xls,ls)-the level structure returned.
!
!  updated parameters-
!  nodlvl-on input, it specifies a section subgraph.
!  on return, it contains the node level numbers.
!
  integer n
!
  integer iadj(*)
  integer j
  integer lbegin
  integer ls(n)
  integer lvl
  integer lvlend
  integer nlvl
  integer node
  integer nodlvl(n)
  integer root
  integer xadj(n+1)
  integer xls(n+1)
!
  call fnroot ( root, xadj, iadj, nodlvl, nlvl, xls, ls, n )

  do lvl = 1, nlvl

    lbegin = xls(lvl)
    lvlend = xls(lvl+1)-1

    do j = lbegin, lvlend
      node = ls(j)
      nodlvl(node) = lvl
    end do

  end do

  return
end
subroutine fnofnz ( xadj, iadj, perm, invprm, nblks, xblk, xnonz, nzsubs, &
  maxnz, n )
!
!*******************************************************************************
!
!! FNOFNZ finds the column subscripts of the off-block-diagonal nonzeros
!  in the lower triangle of a partitioned matrix.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  (perm,invprm)-the permutation vectors.
!  (nblks,xblk)-the block partitioning.
!
!  output parameters-
!  (xnonz,nzsubs)-the column subscripts of the nonzeros
!  of a to the left of the diagonal blocks are
!  stored row by row in contiguous locations in the
!  array nzsubs.  xnonz is the index vector to it.
!
!  updated parameter-
!  maxnz-on input, it contains the size of the vector
!  nzsubs; and on output, the number of nonzeros found.
!
  integer nblks
  integer n
!
  integer iadj(*)
  integer blkbeg
  integer blkend
  integer i
  integer invprm(n)
  integer j
  integer jperm
  integer jxnonz
  integer k
  integer kstop
  integer kstrt
  integer maxnz
  integer nabor
  integer nzcnt
  integer nzsubs(*)
  integer perm(n)
  integer xadj(n+1)
  integer xblk(nblks+1)
  integer xnonz(n+1)
!
  nzcnt = 1

  if ( nblks <= 0 ) then
    maxnz = 0
    return
  end if
!
!  Loop over the blocks.
!
  do i = 1, nblks

    blkbeg = xblk(i)
    blkend = xblk(i+1)-1
!
!  Loop over the rows of the I-th block.
!
    do j = blkbeg, blkend

      xnonz(j) = nzcnt
      jperm = perm(j)
      kstrt = xadj(jperm)
      kstop = xadj(jperm+1)-1
!
!  Loop over the nonzeros of row J.
!
      if ( kstrt <= kstop ) then

        do k = kstrt, kstop

          nabor = iadj(k)
          nabor = invprm(nabor)
!
!  Check to see if it is to the left of the I-th diagonal block.
!
          if ( nabor < blkbeg ) then

            if ( nzcnt <= maxnz ) then
              nzsubs(nzcnt) = nabor
            else
              write ( *, * ) ' '
              write ( *, * ) 'FNOFNZ - Fatal error!'
              write ( *, * ) '  NZCNT exceeds MAXNZ'
              write ( *, * ) '  NZCNT = ', nzcnt
              write ( *, * ) '  MAXNZ = ', maxnz
              stop
            end if

            nzcnt = nzcnt+1

          end if

        end do
!
!  Sort the subscripts of row J.
!
        if ( nzcnt - 1 <= maxnz ) then
          jxnonz = xnonz(j)
          if ( nzcnt > jxnonz ) then
            call sorts1 ( nzcnt-jxnonz, nzsubs(jxnonz) )
          end if
        end if

      end if

    end do

  end do

  xnonz(blkend+1) = nzcnt

  maxnz = nzcnt-1

  return
end
subroutine fnroot ( root, xadj, iadj, mask, nlvl, xls, ls, n )
!
!*******************************************************************************
!
!! FNROOT finds pseudo-peripheral nodes.
!
!
!  Discussion:
!
!    The pseudo-peripheral nodes are find using a modified version of
!    the scheme by gibbs, poole and stockmeyer.  it determines such
!    a node for the section subgraph specified by mask and root.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  mask-specifies a section subgraph. nodes for which
!  mask is zero are ignored by fnroot.
!
!  updated parameter-
!  root-on input, it (along with mask) defines the
!  component for which a pseudo-peripheral node is
!  to be found. on output, it is the node obtained.
!
!  output parameters-
!  nlvl-is the number of levels in the level structure
!  rooted at the node root.
!  (xls,ls)-the level structure array pair containing
!  the level structure found.
!
  integer n
!
  integer iadj(*)
  integer iccsze
  integer j
  integer jstrt
  integer k
  integer kstop
  integer kstrt
  integer ls(n)
  integer mask(n)
  integer mindeg
  integer nabor
  integer ndeg
  integer nlvl
  integer node
  integer nunlvl
  integer root
  integer xadj(n+1)
  integer xls(n+1)
!
!  Determine the level structure rooted at ROOT.
!
  call rootls ( root, xadj, iadj, mask, nlvl, xls, ls, n )

  iccsze = xls(nlvl+1)-1
  if ( nlvl == 1 .or. nlvl == iccsze ) then
    return
  end if
!
!  Pick a node with minimum degree from the last level.
!
   10 continue

  jstrt = xls(nlvl)
  mindeg = iccsze
  root = ls(jstrt)

  if ( iccsze > jstrt ) then

    do j = jstrt, iccsze

      node = ls(j)
      ndeg = 0
      kstrt = xadj(node)
      kstop = xadj(node+1)-1

      do k = kstrt, kstop
        nabor = iadj(k)
        if ( mask(nabor) > 0 ) then
          ndeg = ndeg+1
        end if
      end do

      if ( ndeg < mindeg ) then
        root = node
        mindeg = ndeg
      end if

    end do

  end if
!
!  Generate its rooted level structure.
!
  call rootls ( root, xadj, iadj, mask, nunlvl, xls, ls, n )

  if ( nunlvl <= nlvl ) then
    return
  end if

  nlvl = nunlvl

  if ( nlvl < iccsze ) then
    go to 10
  end if

  return
end
subroutine fnspan ( xadj, iadj, nodlvl, nspan, set, level, nadjs, adjs, &
  leaf, n )
!
!*******************************************************************************
!
!! FNSPAN finds the span of a subset in a level subgraph in a level structure.
!
!
!  Discussion:
!
!    The adjacent set of the span in the lower level is also determined.  
!    If the span has an unnumbered node in the higher level, an unnumbered 
!    leaf node (i.e. one with no neighbor in next level) will be returned.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  level-level number of the current set.
!
!  updated parameters-
!  (nspan,set)-the input set.  on return, it contains
!  the resulting span set.
!  nodlvl-the level number vector.  nodes considered
!  will have their nodlvl changed to zero.
!
!  output parameters-
!  (nadjs,adjs)-the adjacent set of the span in the
!  lower level.
!  leaf-if the span has an unnumbered higher level node,
!  leaf returns an unnumbered leaf node in the level
!  structure, otherwise, leaf is zero.
!
  integer n
!
  integer iadj(*)
  integer adjs(*)
  integer i
  integer j
  integer jstop
  integer jstrt
  integer leaf
  integer level
  integer lvl
  integer lvlm1
  integer nadjs
  integer nbr
  integer nbrlvl
  integer node
  integer nodlvl(n)
  integer nspan
  integer set(*)
  integer setptr
  integer xadj(n+1)
!
!  Initialization.
!
  leaf = 0
  nadjs = 0
  setptr = 0

   10 continue

  setptr = setptr+1

  if ( setptr > nspan ) then
    return
  end if
!
!  For each node in the partially spanned set...
!
  node = set(setptr)
  jstrt = xadj(node)
  jstop = xadj(node+1)-1

  if ( jstop < jstrt ) then
    go to 10
  end if
!
!  For each neighbor of node, test its nodlvl value
!
  do j = jstrt, jstop

    nbr = iadj(j)
    nbrlvl = nodlvl(nbr)

    if ( nbrlvl > 0 ) then

      if ( nbrlvl == level ) then
        go to 30
      end if

      if ( nbrlvl > level ) then
        go to 60
      end if
!
!  NBR is in level-1, add it to ADJS.
!
      nadjs = nadjs+1
      adjs(nadjs) = nbr
      go to 40
!
!  NBR is in level, add it to the span set.
!
   30     continue

      nspan = nspan+1
      set(nspan) = nbr

   40     continue

      nodlvl(nbr) = 0

    end if

  end do

  go to 10
!
!  NBR is in level+1. 
!
!  Find an unnumbered leaf node by tracing a path up the level structure.  
!
!  Then reset the NODLVL values of nodes in adjs.
!
   60 continue

  leaf = nbr
  lvl = level + 1

   70 continue

  jstrt = xadj(leaf)
  jstop = xadj(leaf+1)-1

  do j = jstrt, jstop

    nbr = iadj(j)

    if ( nodlvl(nbr) > lvl ) then
      leaf = nbr
      lvl = lvl + 1
      go to 70
    end if

  end do

  if ( nadjs <= 0 ) then
    return
  end if

  lvlm1 = level-1

  do i = 1, nadjs
    node = adjs(i)
    nodlvl(node) = lvlm1
  end do

  return
end
subroutine fntadj ( xadj, iadj, perm, nblks, xblk, father, n )
!
!*******************************************************************************
!
!! FNTADJ determines the quotient tree adjacency structure for a graph.
!
!
!  Discussion:
!
!    The graph is defined by the matrix that is to be factored.
!
!    The structure is represented by the father vector.
!
!  Modified:
!
!    30 April 2000
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer PERM(N), the permutation vector.
!
!    Input, integer NBLKS, XBLK(NBLKS+1), the tree partitioning.
!
!    Output, integer FATHER(N), the father vector of the quotient tree.
!
!    Input, integer N, the number of equations.
!
  integer nblks
  integer n
!
  integer iadj(*)
  integer father(n)
  integer i
  integer istop
  integer istrt
  integer j
  integer jstop
  integer jstrt
  integer k
  integer mask(n)
  integer nabor
  integer nbrblk
  integer node
  integer perm(n)
  integer xadj(n+1)
  integer xblk(nblks+1)
!
!  Initialize the block number vector.
!
  do k = 1, nblks
    istrt = xblk(k)
    istop = xblk(k+1)-1
    do i = istrt, istop
      node = perm(i)
      mask(node) = k
    end do
  end do
!
!  For each block, find its father block in the tree structure.
!
  father(nblks) = 0

  do k = 1, nblks-1

    istrt = xblk(k)
    istop = xblk(k+1)-1

    do i = istrt, istop

      node = perm(i)
      jstrt = xadj(node)
      jstop = xadj(node+1)-1

      do j = jstrt, jstop

        nabor = iadj(j)
        nbrblk = mask(nabor)

        if ( nbrblk > k ) then
          father(k) = nbrblk
          go to 10
        end if

      end do

    end do

    father(k) = 0

10      continue

  end do

  return
end
subroutine fntenv ( xadj, iadj, perm, invprm, nblks, xblk, xenv, ienv, n )
!
!*******************************************************************************
!
!! FNTENV determines the envelope index vector for the envelope of the
!  diagonal blocks of a tree partitioned system.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer PERM(N), the permutation vector.
!
!    Input, integer INVPRM(N), the inverse permutation vector.
!
!  input parameters-
!
!  (nblks,xblk)-the tree partitioning.
!
!  output parameters-
!  xenv-the envelope index vector.
!  ienv-the size of the envelope found.
!
  integer nblks
  integer n
!
  integer iadj(*)
  integer blkbeg
  integer blkend
  integer i
  integer ienv
  integer ifirst
  integer invprm(n)
  integer j
  integer jstop
  integer jstrt
  integer k
  integer kfirst
  integer nbr
  integer node
  integer perm(n)
  integer xadj(n+1)
  integer xblk(nblks+1)
  integer xenv(n+1)
!
  ienv = 1
!
!  Loop through each block in the partitioning.
!
  do k = 1, nblks

    blkbeg = xblk(k)
    blkend = xblk(k+1)-1
!
!  KFIRST stores the first node in the K-th block that has a neighbor
!  in the previous blocks.
!
    kfirst = blkend

    do i = blkbeg, blkend

      xenv(i) = ienv
      node = perm(i)
      jstrt = xadj(node)
      jstop = xadj(node+1)-1
!
!  IFIRST stores the first nonzero in the I-th row within
!  the K-th block.
!
      ifirst = i

      do j = jstrt, jstop

        nbr = iadj(j)
        nbr = invprm(nbr)

        if ( nbr >= blkbeg ) then
          ifirst = min ( ifirst, nbr )
        else
          ifirst = min ( ifirst, kfirst )
          kfirst = min ( kfirst, i )
        end if

      end do

      ienv = ienv+i-ifirst

    end do

  end do

  xenv(blkend+1) = ienv
  ienv = ienv-1

  return
end
subroutine fn1wd ( root, xadj, iadj, mask, nsep, sep, nlvl, xls, ls, n )
!
!*******************************************************************************
!
!! FN1WD finds one-way dissectors of a connected component.
!
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines (along with MASK) the
!    component to be processed.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameters-
!  nsep-number of nodes in the one-way dissectors.
!  sep-vector containing the dissector nodes.
!
!  updated parameter-
!  mask-nodes in the dissector have their mask values
!  set to zero.
!
!  working parameters-
!  (xls,ls)-level structure used by the routine fnroot.
!
  integer n
!
  integer iadj(*)
  real deltp1
  real fnlvl
  integer i
  integer j
  integer k
  integer kstop
  integer kstrt
  integer lp1beg
  integer lp1end
  integer ls(n)
  integer lvl
  integer lvlbeg
  integer lvlend
  integer mask(n)
  integer nbr
  integer nlvl
  integer node
  integer nsep
  integer root
  integer sep(*)
  real width
  integer xadj(n+1)
  integer xls(n+1)
!
  call fnroot ( root, xadj, iadj, mask, nlvl, xls, ls, n )

  fnlvl = real ( nlvl )
  nsep = xls(nlvl+1)-1
  width = real ( nsep ) / fnlvl
  deltp1 = 1.0 + sqrt ( ( 3.0 * width + 13.0 ) / 2.0 )
!
!  The component is too small, or the level structure
!  is very long and narrow. return the whole component.
!
  if ( nsep < 50 .or. deltp1 > 0.5 * fnlvl ) then

    do i = 1, nsep
      node = ls(i)
      sep(i) = node
      mask(node) = 0
    end do

    return

  end if
!
!  Find the parallel dissectors.
!
  nsep = 0
  i = 0

   30 continue

  i = i + 1
  lvl = int ( real ( i ) * deltp1 + 0.5 )

  if ( lvl >= nlvl ) then
    return
  end if

  lvlbeg = xls(lvl)
  lp1beg = xls(lvl+1)
  lvlend = lp1beg-1
  lp1end = xls(lvl+2)-1

  do j = lp1beg, lp1end
    node = ls(j)
    xadj(node) = -xadj(node)
  end do
!
!  Nodes in level LVL are chosen to form dissector.
!  Include only those with neighbors in lvl+1 level.
!  XADJ is used temporarily to mark nodes in lvl+1.
!
  do j = lvlbeg, lvlend

    node = ls(j)
    kstrt = xadj(node)
    kstop = abs ( xadj(node+1) ) - 1

    do k = kstrt, kstop

      nbr = iadj(k)

      if ( xadj(nbr) <= 0 ) then
        nsep = nsep+1
        sep(nsep) = node
        mask(node) = 0
        go to 60
      end if

    end do

   60   continue

  end do

  do j = lp1beg, lp1end
    node = ls(j)
    xadj(node) = -xadj(node)
  end do

  go to 30
end
subroutine gennd ( n, xadj, iadj, perm, xls, ls )
!
!*******************************************************************************
!
!! GENND finds a nested dissection ordering for a general graph.
!
!
!  Modified:
!
!    30 April 2000
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameters-
!  perm-the nested dissection ordering.
!
!  working parameters-
!  mask-is used to mask off variables that have
!  been numbered during the orderng process.
!  (xls,ls)-this level structure pair is used as
!  temporary storage by fnroot.
!
  integer n
!
  integer iadj(*)
  integer i
  integer ls(n)
  integer mask(n)
  integer nsep
  integer num
  integer perm(n)
  integer root
  integer xadj(n+1)
  integer xls(n+1)
!
!  MASK masks off variables that have been numbered.
!
  mask = 1

  num = 0

  do i = 1, n
!
!  For each masked component...
!
10      continue
!
!  Find a separator and number the nodes next.
!
    if ( mask(i) /= 0 ) then

      root = i

      call fndsep ( root, xadj, iadj, mask, nsep, perm(num+1), xls, ls, n )

      num = num + nsep

      if ( num >= n ) then
        call ivec_reverse ( n, perm )
        return
      end if

      go to 10

    end if

  end do
!
!  Separators found first should be ordered last.
!
  call ivec_reverse ( n, perm )

  return
end
subroutine genqmd ( n, xadj, iadj, perm, invprm, marker, rchset, nbrhd, &
  qsize, qlink, nofsub )
!
!*******************************************************************************
!
!! GENQMD implements the quotient minimum degree algorithm.  
!
!
!  Discussion:
!
!    The routine uses the implicit representation of the elimination 
!    graphs by quotient graphs, and the notion of indistinguishable nodes.
!
!    The adjacency vector IADJ will be destroyed.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameters-
!  perm-the minimum degree ordering.
!  invprm-the inverse of perm.
!
!  working parameters-
!  marker-a marker vector, where marker(i) is
!  negative means node i has been merged with
!  another node and thus can be ignored.
!  rchset-vector used for the reachable set.
!  nbrhd-vector used for the neighborhood set.
!  qsize-vector used to store the size of
!  indistinguishable supernodes.
!  qlink-vector to store indistinguishable nodes,
!  i,qlink(i),qlink(qlink(i)) are the
!  members of the supernode represented by i.
!
  integer n
!
  integer iadj(*)
  integer deg(n)
  integer inode
  integer invprm(n)
  integer ip
  integer irch
  integer j
  integer marker(n)
  integer mindeg
  integer nbrhd(*)
  integer ndeg
  integer nhdsze
  integer node
  integer nofsub
  integer np
  integer num
  integer nump1
  integer nxnode
  integer perm(n)
  integer qlink(n)
  integer qsize(n)
  integer rchset(*)
  integer rchsze
  integer search
  integer thresh
  integer xadj(n+1)
!
!  Initialize degree vector and other working variables.
!
!  DEG(I) negative means node I has been numbered.
!
  mindeg = n
  nofsub = 0

  do node = 1, n
    perm(node) = node
    invprm(node) = node
    marker(node) = 0
    qsize(node) = 1
    qlink(node) = 0
    ndeg = xadj(node+1)-xadj(node)
    deg(node) = ndeg
    mindeg = min ( mindeg, ndeg )
  end do

  num = 0
!
!  Perform threshold search to get a node of minimum degree.
!  Variable search points to where search should start.
!
   20 continue

  search = 1
  thresh = mindeg
  mindeg = n

   30 continue

  nump1 = num+1
  search = max ( search, nump1 )

  do j = search, n

    node = perm(j)

    if ( marker(node) >= 0 ) then
      ndeg = deg(node)
      if ( ndeg <= thresh ) then
        go to 50
      end if
      mindeg = min ( mindeg, ndeg )
    end if

  end do

  go to 20
!
!  The node has minimum degree. 
!  Find its reachable sets by calling qmdrch.
!
   50 continue

  search = j
  nofsub = nofsub+deg(node)
  marker(node) = 1

  call qmdrch ( node, xadj, iadj, deg, marker, rchsze, rchset, nhdsze, nbrhd, n )
!
!  Eliminate all nodes indistinguishable from node.
!  They are given by node, qlink(node),
!
  nxnode = node

   60 continue

  num = num + 1
  np = invprm(nxnode)
  ip = perm(num)
  perm(np) = ip
  invprm(ip) = np
  perm(num) = nxnode
  invprm(nxnode) = num
  deg(nxnode) = -1
  nxnode = qlink(nxnode)

  if ( nxnode > 0 ) then
    go to 60
  end if

  if ( rchsze <= 0 ) then
    go to 80
  end if
!
!  Update the degrees of the nodes in the reachable
!  set and identify indistinguishable nodes.
!
  call qmdupd ( xadj, iadj, rchsze, rchset, deg, qsize, qlink, marker, &
    rchset(rchsze+1), nbrhd(nhdsze+1), n )
!
!  Reset marker value of nodes in reach set.
!  Update threshold value for cyclic search.
!  Call qmdqt to form new quotient graph.
!
  marker(node) = 0

  do irch = 1, rchsze

    inode = rchset(irch)

    if ( marker(inode) >= 0 ) then

      marker(inode) = 0
      ndeg = deg(inode)
      mindeg = min ( mindeg, ndeg )

      if ( ndeg <= thresh ) then
        mindeg = thresh
        thresh = ndeg
        search = invprm(inode)
      end if

    end if

  end do

  if ( nhdsze > 0 ) then
    call qmdqt ( node, xadj, iadj, marker, rchsze, rchset, nbrhd, n )
  end if

   80 continue

  if ( num < n ) then
    go to 30
  end if

  return
end
subroutine genrcm ( n, xadj, iadj, perm, xls )
!
!*******************************************************************************
!
!! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph. 
!
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains 
!    an ordering by calling RCM.
!
!  Modified:
!
!    30 April 2000
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameter-
!  perm-vector that contains the rcm ordering.
!
!  working parameters-
!  xls-the index vector for a level structure.  the
!  level structure is stored in the currently
!  unused spaces in the permutation vector perm.
!
  integer n
!
  integer iadj(*)
  integer i
  integer iccsze
  integer mask(n)
  integer nlvl
  integer num
  integer perm(n)
  integer root
  integer xadj(n+1)
  integer xls(n+1)
!
!  MASK marks variables that have been numbered.
!
  mask = 1

  num = 1

  do i = 1, n
!
!  For each masked connected component...
!
    if ( mask(i) /= 0 ) then

      root = i
!
!  Find a pseudo-peripheral node root.
!  Note that the level structure found by
!  fnroot is stored starting at perm(num).
!
      call fnroot ( root, xadj, iadj, mask, nlvl, xls, perm(num), n )
!
!  RCM orders the component using ROOT as the starting node.
!
      call rcm ( root, xadj, iadj, mask, perm(num), iccsze, xls, n )

      num = num + iccsze

      if ( num > n ) then
        return
      end if

    end if

  end do

  return
end
subroutine genrqt ( n, xadj, iadj, nblks, xblk, perm, xls, ls, nodlvl )
!
!*******************************************************************************
!
!! GENRQT determines a partitioned ordering for a possibly disconnected
!  graph using the refined quotient tree algorithm.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameters-
!  (nblks,xblk)-the quotient tree partitioning.
!  perm-the permutation vector.
!
!  working parameters-
!  (xls,ls)-this level structure pair is used by
!  fnroot to find a pseudo-peripheral node.
!  nodlvl-a temporary vector to store the level
!  number of each node in a level structure.
!
  integer n
!
  integer iadj(*)
  integer i
  integer ixls
  integer leaf
  integer ls(n)
  integer nblks
  integer nlvl
  integer nodlvl(n)
  integer perm(n)
  integer root
  integer xadj(n+1)
  integer xblk(*)
  integer xls(n+1)
!
  do i = 1, n
    nodlvl(i) = 1
  end do

  nblks = 0
  xblk(1) = 1

  do i = 1, n
!
!  For each connected component...
!
    if ( nodlvl(i) > 0 ) then
!
!  Find a rooted level structure.
!
      root = i
      call fnlvls ( root, xadj, iadj, nodlvl, nlvl, xls, ls, n )
      ixls = xls(nlvl)
      leaf = ls(ixls)
!
!  RQTREE gets the block order.
!
      call rqtree ( leaf, xadj, iadj, perm, nblks, xblk, nodlvl, xls, ls, n )

    end if

  end do

  return
end
subroutine gen1wd ( n, xadj, iadj, nblks, xblk, perm, xls, ls )
!
!*******************************************************************************
!
!! GEN1WD finds a one-way dissection partitioning for a general graph.
!
!
!  Discussion:
!
!    FN1WD is called once for each connected component.
!
!  Modified:
!
!    30 April 2000
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the number of equations.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer NBLKS, integer XBLK(*), the partitioning found.
!
!    Output, integer PERM(N), the one-way dissection ordering.
!
!    Workspace, integer XLS(N+1), LS(N), a level structure used
!    by ROOTLS.
!
  integer n
!
  integer iadj(*)
  integer i
  integer iccsze
  integer j
  integer k
  integer lnum
  integer ls(n)
  integer mask(n)
  integer nblks
  integer nlvl
  integer node
  integer nsep
  integer num
  integer perm(n)
  integer root
  integer xadj(n+1)
  integer xblk(*)
  integer xls(n+1)
!
!  MASK will mark variables that have been numbered.
!
  mask = 1

  nblks = 0
  num = 0

  do i = 1, n

    if ( mask(i) /= 0 ) then
!
!  Find a one-way dissector for each component.
!
      root = i

      call fn1wd ( root, xadj, iadj, mask, nsep, perm(num+1), nlvl, xls, ls, n )

      num = num + nsep
      nblks = nblks + 1
      xblk(nblks) = n - num + 1
      iccsze = xls(nlvl+1)-1
!
!  Number the remaining nodes in the component.
!  Each component in the remaining subgraph forms
!  a new block in the partitioning.
!
      do j = 1, iccsze

        node = ls(j)

        if ( mask(node) /= 0 ) then

          call rootls ( node, xadj, iadj, mask, nlvl, xls, perm(num+1), n )

          lnum = num + 1
          num = num + xls(nlvl+1) - 1
          nblks = nblks+1
          xblk(nblks) = n - num + 1

          do k = lnum, num
            node = perm(k)
            mask(node) = 0
          end do

          if ( num > n ) then
            go to 50
          end if

        end if

      end do

    end if

  end do
!
!  Since dissectors found first should be ordered last,
!  REVRSE is called to adjust the ordering
!  vector and the block index vector.
!
   50 continue

  call ivec_reverse ( n, perm )

  call ivec_reverse ( nblks, xblk )

  xblk(nblks+1) = n + 1

  return
end
subroutine gs_factor ( n, ixlnz, lnz, xnzsub, nzsub, diag, link, first )
!
!*******************************************************************************
!
!! GS_FACTOR performs the symmetric factorization for a general sparse
!  system, stored in the compressed subscript data format.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!  ixlnz-index vector for lnz.  ixlnz(i) points to the
!  start of nonzeros in column i of factor l.
!  (xnzsub,nzsub)-the compressed subscript data
!  structure for factor l.
!
!  updated parameters-
!  lnz-on input,contains nonzeros of a,and on
!  return,the nonzeros of l.
!  diag-the diagonal of l overwrites that of a.
!
!  working parameters-
!  link-at step j,the list in
!  link(j),link(link(j)),.
!  consists of those columns that will modify
!  the column l(*,j).
!  first-temporary vector to point to the first
!  nonzero in each column that will be used
!  next for modification.
!
  integer n
!
  real diag(n)
  real diagj
  integer first(n)
  integer i
  integer ii
  integer istop
  integer istrt
  integer isub
  integer ixlnz(n+1)
  integer j
  integer k
  integer kfirst
  integer link(n)
  real ljk
  real lnz(*)
  integer newk
  integer nzsub(*)
  real temp(n)
  integer xnzsub(*)
!
!  Initialize working vectors.
!
  link = 0
  temp = 0.0
!
!  Compute column l(*,j) for j = 1,....,N.
!
  do j = 1, n
!
!  For each column l(*,k) that affects l(*,j).
!
    diagj = 0.0
    newk = link(j)

   20   continue

    k = newk
    if ( k == 0 ) then
      go to 40
    end if
    newk = link(k)
!
!  Outer product modification of L(*,J) by L(*,K) starting at FIRST(K) of L(*,K).
!
    kfirst = first(k)
    ljk = lnz(kfirst)
    diagj = diagj + ljk*ljk
    istrt = kfirst+1
    istop = ixlnz(k+1)-1

    if ( istop < istrt ) then
      go to 20
    end if
!
!  Before modification, update vectors first,
!  and link for future modification steps.
!
    first(k) = istrt
    i = xnzsub(k) + ( kfirst - ixlnz(k) ) + 1
    isub = nzsub(i)
    link(k) = link(isub)
    link(isub) = k
!
!  The actual mod is saved in vector temp.
!
    do ii = istrt, istop
      isub = nzsub(i)
      temp(isub) = temp(isub) + lnz(ii) * ljk
      i = i + 1
    end do

    go to 20
!
!  Apply the modifications accumulated in temp to column L(*,J).
!
   40   continue

    diagj = diag(j) - diagj

    if ( diagj <= 0.0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'GSFCT - Fatal error'
      write ( *, * ) '  Zero or negative diagonal entry!'
      write ( *, * ) '  DIAG(J) = ', diagj
      write ( *, * ) '  for diagonal J = ', j
      write ( *, * ) ' '
      write ( *, * ) '  The matrix is not positive definite.'
      stop
    end if

    diagj = sqrt ( diagj )
    diag(j) = diagj
    istrt = ixlnz(j)
    istop = ixlnz(j+1) - 1

    if ( istop >= istrt ) then

      first(j) = istrt
      i = xnzsub(j)
      isub = nzsub(i)
      link(j) = link(isub)
      link(isub) = j

      do ii = istrt, istop
        isub = nzsub(i)
        lnz(ii) = ( lnz(ii) - temp(isub) ) / diagj
        temp(isub) = 0.0
        i = i+1
      end do

    end if

  end do

  return
end
subroutine gs_solve ( n, ixlnz, lnz, xnzsub, nzsub, diag, rhs )
!
!*******************************************************************************
!
!! GS_SOLVE solves a factored system, stored in compressed subscript format.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!  (ixlnz,lnz)-structure of nonzeros in l.
!  (xnzsub,nzsub)-compressed subscript structure.
!  diag-diagonal components of l.
!
!  updated parameter-
!  rhs-on input,it contains the rhs vector,and on
!  output, the solution vector.
!
  integer n
!
  real diag(n)
  integer i
  integer ii
  integer istop
  integer istrt
  integer isub
  integer ixlnz(n+1)
  integer j
  integer jj
  real lnz(*)
  integer nzsub(*)
  real rhs(n)
  real rhsj
  real s
  integer xnzsub(*)
!
!  Forward substitution.
!
  do j = 1, n

    rhsj = rhs(j) / diag(j)
    rhs(j) = rhsj
    istrt = ixlnz(j)
    istop = ixlnz(j+1) - 1

    i = xnzsub(j)

    do ii = istrt, istop
      isub = nzsub(i)
      rhs(isub) = rhs(isub) - lnz(ii) * rhsj
      i = i + 1
    end do

  end do
!
!  Backward substitution.
!
  j = n

  do jj = 1, n

    s = rhs(j)
    istrt = ixlnz(j)
    istop = ixlnz(j+1) - 1

    i = xnzsub(j)

    do ii = istrt, istop
      isub = nzsub(i)
      s = s - lnz(ii) * rhs(isub)
      i = i + 1
    end do

    rhs(j) = s / diag(j)
    j = j - 1

  end do

  return
end
subroutine i_swap ( i, j )
!
!*******************************************************************************
!
!! I_SWAP switches two integer values.
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
subroutine ivec_copy ( n, a, b )
!
!*******************************************************************************
!
!! IVEC_COPY copies one integer vector into another.
!
!
!  Parameters:
!
!    Input, integer N, the number of entries in A.
!
!    Input, integer A(N), the vector to be copied.
!
!    Output, integer B(N), a copy of A.
!
  integer n
!
  integer a(n)
  integer b(n)
!
  b = a

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
subroutine perm_inverse ( n, perm, invprm )
!
!*******************************************************************************
!
!! PERM_INVERSE produces the inverse of a given permutation.
!
!
!  Modified:
!
!    29 April 2000
!
!  Parameters:
!
!    Input, integer N, the number of equations.
!
!    Input, integer PERM(N), contains the reordering of the
!    variables and equations.
!
!    Output, integer INVPRM(N), the inverse ordering with the
!    property that INVPRM(PERM(I)) = I.
!
  integer n
!
  integer i
  integer invprm(n)
  integer k
  integer perm(n)
!
  do i = 1, n
    k = perm(i)
    invprm(k) = i
  end do

  return
end
subroutine perm_rv ( n, rhs, perm )
!
!*******************************************************************************
!
!! PERM_RV should be called once the linear system has been solved and
!  the solution returned in RHS.  The routine then undoes the permutation
!  of RHS, restoring the original ordering.  To do this, it needs the
!  PERM vector which defined the reordering used by the solver.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the number of equations.
!
!    Input/output, real RHS(N).
!    On input, RHS contains the solution of the permuted linear system.
!    On output, RHS contains the solution of the original linear system.
!
!    Input, integer PERM(N), the permutation information.
!    PERM(I) = K means that the K-th equation and variable in the
!    original ordering became the I-th equation and variable in the
!    reordering.
!
  integer n
!
  integer i
  integer iput
  integer istart
  integer perm(n)
  real pull
  real put
  real rhs(n)
!
!  Mark PERM with negative signs which will be removed
!  as each permuted element is restored to its rightful place
!
  do i = 1, n
    perm(i) = -perm(i)
  end do
!
!  Search for the next element of perm which is the first
!  element of a permutation cycle
!
  istart = 0

   20 continue

  istart = istart+1

  if ( istart > n ) then
    return
  end if

  if ( perm(istart) > 0 ) then
    go to 20
  end if

  if ( abs ( perm(istart) ) /= istart ) then
    go to 30
  end if

  perm(istart) = abs ( perm(istart) )
  go to 20
!
!  Begin a cycle.
!
   30 continue

  perm(istart) = abs ( perm(istart) )
  iput = istart
  pull = rhs(iput)

   40 continue

  iput = abs ( perm(iput) )
  put = rhs(iput)
  rhs(iput) = pull
  pull = put

  if ( perm(iput) > 0 ) then
    go to 20
  end if

  perm(iput) = abs ( perm(iput) )
  go to 40

end
subroutine qmdmrg ( xadj, iadj, deg, qsize, qlink, marker, deg0, nhdsze, &
  nbrhd, rchset, ovrlp, n )
!
!*******************************************************************************
!
!! QMDMRG merges indistinguishable nodes in the minimum degree ordering
!  algorithm, and computes the new degrees of these new supernodes.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  deg0-the number of nodes in the given set.
!  (nhdsze,nbrhd)-the set of eliminated supernodes
!  adjacent to some nodes in the set.
!
!  updated parameters-
!  deg-the degree vector.
!  qsize-size of indistinguishable nodes.
!  qlink-linked list for indistinguishable nodes.
!  marker-the given set is given by those nodes with
!  marker value set to 1.  those nodes with degree
!  updated will have marker value set to 2.
!
!  working parameters-
!  rchset-the reachable set.
!  ovrlp- temp vector to store the intersection of two
!  reachable sets.
!
  integer n
!
  integer iadj(*)
  integer deg(n)
  integer deg0
  integer deg1
  integer head
  integer inhd
  integer iov
  integer irch
  integer j
  integer jstrt
  integer jstop
  integer link
  integer lnode
  integer mark
  integer marker(n)
  integer mrgsze
  integer nabor
  integer nbrhd(*)
  integer nhdsze
  integer node
  integer novrlp
  integer ovrlp(*)
  integer qlink(n)
  integer qsize(n)
  integer rchset(*)
  integer rchsze
  integer root
  integer xadj(n+1)
!
  if ( nhdsze <= 0 ) then
    return
  end if

  do inhd = 1, nhdsze
    root = nbrhd(inhd)
    marker(root) = 0
  end do
!
!  Loop through each eliminated supernode in the set (nhdsze,nbrhd).
!
  do inhd = 1, nhdsze

    root = nbrhd(inhd)
    marker(root) = - 1
    rchsze = 0
    novrlp = 0
    deg1 = 0

   20   continue

    jstrt = xadj(root)
    jstop = xadj(root+1)-1
!
!  Determine the reachable set and its intersection with the input
!  reachable set.
!
    do j = jstrt, jstop

      nabor = iadj(j)
      root = -nabor

      if ( nabor < 0 ) then
        go to 20
      end if

      if ( nabor == 0 ) then
        go to 70
      end if

      mark = marker(nabor)

      if ( mark == 0 ) then

        rchsze = rchsze+1
        rchset(rchsze) = nabor
        deg1 = deg1+qsize(nabor)
        marker(nabor) = 1

      else if ( mark == 1 ) then

        novrlp = novrlp+1
        ovrlp(novrlp) = nabor
        marker(nabor) = 2

      end if

    end do
!
!  From the overlapped set, determine the nodes that can be merged.
!
   70   continue

    head = 0
    mrgsze = 0

    do iov = 1, novrlp

      node = ovrlp(iov)
      jstrt = xadj(node)
      jstop = xadj(node+1)-1

      do j = jstrt, jstop

        nabor = iadj(j)

        if ( marker(nabor) == 0 ) then
          marker(node) = 1
          go to 110
        end if

      end do
!
!  NODE belongs to the new merged supernode.
!  Update the vectors QLINK and QSIZE.
!
      mrgsze = mrgsze + qsize(node)
      marker(node) = -1
      lnode = node

   90     continue

      link = qlink(lnode)

      if ( link > 0 ) then
        lnode = link
        go to 90
      end if

      qlink(lnode) = head
      head = node

  110     continue

    end do

    if ( head > 0 ) then
      qsize(head) = mrgsze
      deg(head) = deg0+deg1-1
      marker(head) = 2
    end if
!
!  Reset marker values.
!
    root = nbrhd(inhd)
    marker(root) = 0

    do irch = 1, rchsze
      node = rchset(irch)
      marker(node) = 0
    end do

  end do

  return
end
subroutine qmdqt ( root, xadj, iadj, marker, rchsze, rchset, nbrhd, n )
!
!*******************************************************************************
!
!! QMDQT performs the quotient graph transformation after a node has
!  been eliminated.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node just eliminated.  It becomes the
!    representative of the new supernode.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  (rchsze,rchset)-the reachable set of root in the
!  old quotient graph.
!  nbrhd-the neighborhood set which will be merged
!  with root to form the new supernode.
!  marker-the marker vector.
!
!  updated parameter-
!  iadj-becomes the adjacency structure of the quotient graph.
!
  integer n
!
  integer iadj(*)
  integer inhd
  integer irch
  integer j
  integer jstrt
  integer jstop
  integer link
  integer marker(n)
  integer nabor
  integer nbrhd(*)
  integer node
  integer rchset(*)
  integer rchsze
  integer root
  integer xadj(n+1)
!
  irch = 0
  inhd = 0
  node = root

   10 continue

  jstrt = xadj(node)
  jstop = xadj(node+1)-2
!
!  Place reach nodes into the adjacent list of node.
!
  do j = jstrt, jstop
    irch = irch+1
    iadj(j) = rchset(irch)
    if ( irch >= rchsze ) then
      go to 40
    end if
  end do
!
!  Link to other space provided by the nbrhd set.
!
  link = iadj(jstop+1)
  node = -link

  if ( link >= 0 ) then
    inhd = inhd + 1
    node = nbrhd(inhd)
    iadj(jstop+1) = -node
  end if

  go to 10
!
!  All reachable nodes have been saved.  End the adj list.
!  Add root to the nbr list of each node in the reach set.
!
   40 continue

  iadj(j+1) = 0

  do irch = 1, rchsze

    node = rchset(irch)

    if ( marker(node) >= 0 ) then

      jstrt = xadj(node)
      jstop = xadj(node+1)-1

      do j = jstrt, jstop

        nabor = iadj(j)

        if ( marker(nabor) < 0 ) then
          iadj(j) = root
          go to 60
        end if

      end do

    end if

60      continue

  end do

  return
end
subroutine qmdrch ( root, xadj, iadj, deg, marker, rchsze, rchset, nhdsze, &
  nbrhd, n )
!
!*******************************************************************************
!
!! QMDRCH determines the reachable set of a node through a given subset.  
!
!
!  Discussion:
!
!    The adjacency structure is assumed to be stored in a quotient graph format.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the given node, which is not in the subset.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  deg-the degree vector.  deg(i) lt 0 means the node
!  belongs to the given subset.
!
!  output parameters-
!  (rchsze,rchset)-the reachable set.
!  (nhdsze,nbrhd)-the neighborhood set.
!
!  updated parameters-
!  marker-the marker vector for reach and nbrhd sets.
!  gt 0 means the node is in reach set.
!  lt 0 means the node has been merged with
!  others in the quotient or it is in nbrhd set.
!
  integer n
!
  integer iadj(*)
  integer deg(n)
  integer i
  integer istop
  integer istrt
  integer j
  integer jstop
  integer jstrt
  integer marker(n)
  integer nabor
  integer nbrhd(*)
  integer nhdsze
  integer node
  integer rchset(*)
  integer rchsze
  integer root
  integer xadj(n+1)
!
!  Loop through the neighbors of ROOT in the quotient graph.
!
  nhdsze = 0
  rchsze = 0
  istrt = xadj(root)
  istop = xadj(root+1)-1

  do i = istrt, istop

    nabor = iadj(i)
    if ( nabor == 0 ) then
      return
    end if

    if ( marker(nabor) /= 0 ) then
      go to 50
    end if
!
!  Include NABOR in the reachable set.
!
    if ( deg(nabor) >= 0 ) then

      rchsze = rchsze+1
      rchset(rchsze) = nabor
      marker(nabor) = 1
      go to 50

    end if
!
!  NABOR has been eliminated.  Find nodes reachable from it.
!
    marker(nabor) = -1
    nhdsze = nhdsze+1
    nbrhd(nhdsze) = nabor

   20   continue

    jstrt = xadj(nabor)
    jstop = xadj(nabor+1)-1

    do j = jstrt, jstop

      node = iadj(j)
      nabor = -node

      if ( node < 0 ) then
        go to 20
      end if

      if ( node == 0 ) then
        go to 50
      end if

      if ( marker(node) == 0 ) then
        rchsze = rchsze+1
        rchset(rchsze) = node
        marker(node) = 1
      end if

    end do

   50   continue

  end do

  return
end
subroutine qmdupd ( xadj, iadj, nlist, list, deg, qsize, qlink, marker, &
  rchset, nbrhd, n )
!
!*******************************************************************************
!
!! QMDUPD performs degree update for a set of nodes in the minimum
!  degree algorithm.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  (nlist,list)-the list of nodes whose degree has to
!  be updated.
!
!  updated parameters-
!  deg-the degree vector.
!  qsize-size of indistinguishable supernodes.
!  qlink-linked list for indistinguishable nodes.
!  marker-used to mark those nodes in reach/nbrhd sets.
!
!  working parameters-
!  rchset-the reachable set.
!  nbrhd- the neighborhood set.
!
  integer n
!
  integer iadj(*)
  integer deg(n)
  integer deg0
  integer deg1
  integer il
  integer inhd
  integer inode
  integer irch
  integer j
  integer jstop
  integer jstrt
  integer list(*)
  integer mark
  integer marker(n)
  integer nabor
  integer nbrhd(*)
  integer nhdsze
  integer nlist
  integer node
  integer qlink(n)
  integer qsize(n)
  integer rchset(*)
  integer rchsze
  integer xadj(n+1)
!
!  Find all eliminated supernodes that are adjacent to some nodes in the
!  given list.  Put them into (nhdsze,nbrhd).  DEG0 contains the number of
!  nodes in the list.
!
  if ( nlist <= 0 ) then
    return
  end if

  deg0 = 0
  nhdsze = 0

  do il = 1, nlist

    node = list(il)
    deg0 = deg0+qsize(node)
    jstrt = xadj(node)
    jstop = xadj(node+1)-1

    do j = jstrt, jstop

      nabor = iadj(j)

      if ( marker(nabor) == 0 .and. deg(nabor) < 0 ) then
        marker(nabor) = -1
        nhdsze = nhdsze+1
        nbrhd(nhdsze) = nabor
      end if

    end do

  end do
!
!  Merge indistinguishable nodes in the list by calling QMDMRG.
!
  if ( nhdsze > 0 ) then
    call qmdmrg ( xadj, iadj, deg, qsize, qlink, marker, deg0, nhdsze, &
      nbrhd, rchset, nbrhd(nhdsze+1), n )
  end if
!
!  Find the new degrees of the nodes that have not been merged.
!
  do il = 1, nlist

    node = list(il)
    mark = marker(node)

    if ( mark == 0 .or. mark == 1 ) then

      marker(node) = 2

      call qmdrch ( node, xadj, iadj, deg, marker, rchsze, rchset, nhdsze, &
        nbrhd, n )

      deg1 = deg0

      do irch = 1, rchsze
        inode = rchset(irch)
        deg1 = deg1 + qsize(inode)
        marker(inode) = 0
      end do

      deg(node) = deg1-1

      do inhd = 1, nhdsze
        inode = nbrhd(inhd)
        marker(inode) = 0
      end do

    end if

  end do

  return
end
subroutine rcm ( root, xadj, iadj, mask, perm, iccsze, deg, n )
!
!*******************************************************************************
!
!! RCM numbers a connected component using the reverse Cuthill McKee algorithm.
!
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected component.
!    It is used as the starting point for the RCM ordering.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  updated parameters-
!  mask-only those nodes with nonzero input mask
!  values are considered by the routine.  the
!  nodes numbered by rcm will have their
!  mask values set to zero.
!
!  output parameters-
!  perm-will contain the rcm ordering.
!  iccsze-is the size of the connected component
!  that has been numbered by rcm.
!
!  working parameter-
!  deg-is a temporary vector used to hold the degree
!  of the nodes in the section graph specified
!  by mask and root.
!
  integer n
!
  integer iadj(*)
  integer deg(n)
  integer fnbr
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer k
  integer l
  integer lbegin
  integer lnbr
  integer lperm
  integer lvlend
  integer mask(n)
  integer nbr
  integer node
  integer perm(n)
  integer root
  integer xadj(n+1)
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
  call degree ( root, xadj, iadj, mask, deg, iccsze, perm, n )

  mask(root) = 0

  if ( iccsze <= 1 ) then
    return
  end if

  lvlend = 0
  lnbr = 1
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
   10 continue

  lbegin = lvlend + 1
  lvlend = lnbr

  do i = lbegin, lvlend
!
!  For each node in current level...
!
    node = perm(i)
    jstrt = xadj(node)
    jstop = xadj(node+1)-1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last unnumbered neighbors 
!  of the current node in PERM.
!
    fnbr = lnbr+1

    do j = jstrt, jstop

      nbr = iadj(j)

      if ( mask(nbr) /= 0 ) then
        lnbr = lnbr+1
        mask(nbr) = 0
        perm(lnbr) = nbr
      end if

    end do

    if ( fnbr >= lnbr ) then
      go to 60
    end if
!
!  Sort the neighbors of node in increasing order by degree.
!  Linear insertion is used.
!
    k = fnbr

   30   continue

    l = k
    k = k+1
    nbr = perm(k)

   40   continue

    if ( l > fnbr ) then

      lperm = perm(l)

      if ( deg(lperm) > deg(nbr) ) then
        perm(l+1) = lperm
        l = l-1
        go to 40
      end if

    end if

    perm(l+1) = nbr
    if ( k < lnbr ) then
      go to 30
    end if

   60   continue

  end do

  if ( lnbr > lvlend ) then
    go to 10
  end if
!
!  We now have the Cuthill-McKee ordering.  Reverse it.
!
  call ivec_reverse ( iccsze, perm )

  return
end
subroutine reach ( root, xadj, iadj, mask, marker, rchsze, rchset, nhdsze, &
  nbrhd, n )
!
!*******************************************************************************
!
!! REACH determines the reachable set of a node through a subset in a subgraph.
!
!
!  Discussion:
!
!    The routine returns the neighborhood set of node Y in subset S, that is,
!    NBRHD(Y,S), the set of nodes in S that can be reached from Y through
!    a subset of S.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the given node, which is not in the subset S.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  mask-the mask vector for the set s.
!   = 0, if the node is not in s,
!  > 0, if the node is in s.
!
!  output parameters-
!  (nhdsze,nbrhd)-the neighborhood set.
!  (rchsze,rchset)-the reachable set.
!
!  updated parameters-
!  marker-the marker vector used to define the subgraph,
!  nodes in the subgraph have marker value 0.
!  on return, the reachable and neighborhood node
!  sets have their marker values reset to root.
!
  integer n
!
  integer iadj(*)
  integer i
  integer istop
  integer istrt
  integer j
  integer jstop
  integer jstrt
  integer marker(n)
  integer mask(n)
  integer nabor
  integer nbr
  integer nbrhd(*)
  integer nhdptr
  integer nhdsze
  integer node
  integer rchset(*)
  integer rchsze
  integer root
  integer xadj(n+1)
!
!  Initialization.
!
  nhdsze = 0
  rchsze = 0

  if ( marker(root) <= 0 ) then
    rchsze = 1
    rchset(1) = root
    marker(root) = root
  end if

  istrt = xadj(root)
  istop = xadj(root+1)-1
  if ( istop < istrt ) then
    return
  end if
!
!  Loop through the neighbors of ROOT.
!
  do i = istrt, istop

    nabor = iadj(i)

    if ( marker(nabor) /= 0 ) then
      go to 60
    end if
!
!  If NABOR is not in subset S, include it in the reach set.
!
    if ( mask(nabor) <=  0 ) then
      rchsze = rchsze + 1
      rchset(rchsze) = nabor
      marker(nabor) = root
      go to 60
    end if
!
!  NABOR is in subset S, and has not been considered.
!  Include it into the neighborhood set and find the nodes
!  reachable from ROOT through NABOR.
!
   20   continue

    nhdsze = nhdsze+1
    nbrhd(nhdsze) = nabor
    marker(nabor) = root
    nhdptr = nhdsze

   30   continue

    node = nbrhd(nhdptr)
    jstrt = xadj(node)
    jstop = xadj(node+1)-1

    do j = jstrt, jstop

      nbr = iadj(j)

      if ( marker(nbr) == 0 ) then

        if ( mask(nbr) /= 0 ) then
          nhdsze = nhdsze + 1
          nbrhd(nhdsze) = nbr
          marker(nbr) = root
        else
          rchsze = rchsze+1
          rchset(rchsze) = nbr
          marker(nbr) = root
        end if

      end if

    end do

    nhdptr = nhdptr + 1

    if ( nhdptr <= nhdsze ) then
      go to 30
    end if

   60   continue

  end do

  return
end
subroutine rootls ( root, xadj, iadj, mask, nlvl, xls, ls, n )
!
!*******************************************************************************
!
!! ROOTLS generates the connected level structure rooted at a given node.
!
!
!  Discussion:
!
!    Only nodes for which MASK is nonzero will be considered.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node at which the level structure
!    is to be rooted.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(N).
!
!    On input, MASK is 0 for those nodes which should be ignored.
!    On output, every node which was included in the new level
!    structure has had its entry of MASK set to zero.
!
!    Output, integer NLVL, the number of levels in the level
!    structure.  ROOT is in level 1.  The neighbors of ROOT
!    are in level 2, and so on.
!
!    Output, integer XLS(N+1), LS(*).
!    The array pair for the rooted level structure.
!
!    Input, integer N, the number of equations.
!
  integer n
!
  integer iadj(*)
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer ls(n)
  integer lvlend
  integer lvsize
  integer mask(n)
  integer nbr
  integer nlvl
  integer node
  integer root
  integer xadj(n+1)
  integer xls(n+1)
!
  mask(root) = 0
  ls(1) = root
  nlvl = 0
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
   10 continue

  lbegin = lvlend + 1
  lvlend = iccsze
  nlvl = nlvl + 1
  xls(nlvl) = lbegin
!
!  Generate the next level by finding all the masked neighbors of nodes
!  in the current level.
!
  do i = lbegin, lvlend

    node = ls(i)
    jstrt = xadj(node)
    jstop = xadj(node+1)-1

    do j = jstrt, jstop

      nbr = iadj(j)

      if ( mask(nbr) /= 0 ) then
        iccsze = iccsze + 1
        ls(iccsze) = nbr
        mask(nbr) = 0
      end if

    end do

  end do
!
!  Compute the current level width.
!  If it is nonzero, generate the next level.
!
  lvsize = iccsze-lvlend

  if ( lvsize > 0 ) then
    go to 10
  end if
!
!  Reset MASK to one for the nodes in the level structure.
!
  xls(nlvl+1) = lvlend + 1

  do i = 1, iccsze
    node = ls(i)
    mask(node) = 1
  end do

  return
end
subroutine rqtree ( leaf, xadj, iadj, perm, nblks, xblk, nodlvl, adjs, &
  stack, n )
!
!*******************************************************************************
!
!! RQTREE finds a quotient tree ordering for the component specified
!  by leaf and nodlvl.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer LEAF, the input node that defines the connected
!    component.  It is also a leaf node in the rooted level structure 
!    passed to RQTREE, that is, it has no neighbor in the next level.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  output parameters-
!  perm-the permutation vector containing the ordering.
!  (nblks,xblk)-the quotient tree partitioning.
!
!  updated parameters-
!  nodlvl-the node level number vector.  nodes in the
!  component have their nodlvl set to zero as
!  they are numbered.
!
!  working parameters-
!  adjs-temporary vector to store the adjacent set
!  of nodes in a particular level.
!  stack-temporary vector used to maintain the stack
!  of node subsets.  it is organised as-
!  (subset nodes,subset size,subset level).
!
  integer nblks
  integer n
!
  integer iadj(*)
  integer adjs(*)
  integer blksze
  integer ip
  integer j
  integer jp
  integer leaf
  integer level
  integer nadjs
  integer node
  integer nodlvl(n)
  integer npop
  integer nuleaf
  integer num
  integer perm(n)
  integer stack(*)
  integer toplvl
  integer topstk
  integer xadj(n+1)
  integer xblk(*)
!
!  Initialize the stack vector and its pointers.
!
  stack(1) = 0
  stack(2) = 0
  topstk = 2
  toplvl = 0
  num = xblk(nblks+1)-1
!
!  Form a leaf block, that is, one with no neighbors
!  in its next higher level.
!
10    continue

  level = nodlvl(leaf)
  nodlvl(leaf) = 0
  perm(num+1) = leaf
  blksze = 1

  call fnspan ( xadj, iadj, nodlvl, blksze, perm(num+1), level, nadjs, adjs, &
    nuleaf, n )

  if ( nuleaf <= 0 ) then
    go to 30
  end if

  jp = num

  do j = 1, blksze
    jp = jp+1
    node = perm(jp)
    nodlvl(node) = level
  end do

  leaf = nuleaf
  go to 10
!
!  A new block has been found.
!
   30 continue

  nblks = nblks+1
  xblk(nblks) = num+1
  num = num+blksze
!
!  Find the next possible block by using the adjacent
!  set in the lower level and the top node subset (if
!  appropriate) in the stack.
!
  level = level-1

  if ( level <= 0 ) then
    go to 50
  end if

  call ivec_copy ( nadjs, adjs, perm(num+1) )
  blksze = nadjs
!
!  The level of the node subset at the top of the
!  stack is the same as that of the adjacent set.
!  Pop the node subset from the stack.
!
  if ( level == toplvl ) then
    npop = stack(topstk-1)
    topstk = topstk-npop-2
    ip = num + blksze + 1
    call ivec_copy ( npop, stack(topstk+1), perm(ip) )
    blksze = blksze + npop
    toplvl = stack(topstk)
  end if

  call fnspan ( xadj, iadj, nodlvl, blksze, perm(num+1), level, nadjs, &
    adjs, nuleaf, n )

  if ( nuleaf <= 0 ) then
    go to 30
  end if
!
!  Push the current node set into the stack.
!
  call ivec_copy ( blksze, perm(num+1), stack(topstk+1) )
  topstk = topstk+blksze+2
  stack(topstk-1) = blksze
  stack(topstk) = level
  toplvl = level
  leaf = nuleaf
  go to 10
!
!  Before exit.
!
   50 continue

  xblk(nblks+1) = num + 1

  return
end
subroutine shomat ( iadj, iband, ienv, invprm, nadj, n, perm, xadj )
!
!*******************************************************************************
!
!! SHOMAT displays a symbolic picture of a matrix.
!
!
!  Discussion:
!
!    The matrix is defined by the adjacency information in xadj and iadj, 
!    with a possible permutation through perm and invprm.  The routine 
!    also computes the bandwidth and the size of the envelope.
!
!    If no permutation has been done, you must set invprm(i) = perm(i) = i
!    before calling shomat.  otherwise, you must call PERM_INVERSE to
!    get the inverse permutation invprm before calling showmat.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer IBAND, the bandwidth of the matrix.
!
!  ienv   - computed by shomat, the number of cells in the envelope.
!           you could think of this number as the sum of the
!           bandwidths of each row.
!
!    Input, integer INVPRM(N), the inverse permutation.
!
!    Input, integer NADJ, the number of adjacency entries in IADJ.
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
  integer nadj
  integer n
!
  integer iadj(nadj)
  character*1 band(100)
  integer i
  integer iadd
  integer iband
  integer ienv
  integer invprm(n)
  integer itemp
  integer j
  integer jhi
  integer jlo
  integer k
  integer nonz
  integer perm(n)
  integer xadj(n+1)
!
  iband = 0
  ienv = 0
  nonz = 0

  if ( n > 100 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SHOMAT - Fatal error!'
    write ( *, * ) '  N is too large!'
    write ( *, * ) '  Maximum legal value is 100.'
    write ( *, * ) '  Your input value was ', n
    stop
  end if

  write ( *, * ) ' '
  write ( *, * ) 'SHOMAT - Display nonzero structure of matrix.'
  write ( *, * ) ' '

    do i = 1, n

      do k = 1, n
        band(k) = ' '
      end do

      band(i) = 'X'

      iadd = 0
      jlo = xadj(perm(i))
      jhi = xadj(perm(i)+1)-1

      do j = jlo, jhi
        itemp = invprm(iadj(j))
        if ( itemp < i ) then
          nonz = nonz + 1
        end if
        iband = max ( iband, i-itemp )
        band(itemp) = 'X'
        iadd = max ( iadd, i-itemp )
      end do

      write ( *, '(i6,1x,100a1)' ) i, ( band(j), j = 1, n )
      ienv = ienv + iadd

    end do

  write ( *, * ) ' '
  write ( *, * ) 'Lower bandwidth = ', iband
  write ( *, * ) 'Lower envelope contains ', nonz, ' nonzeros.'
  write ( *, * ) 'Lower envelope contains ', ienv, ' entries.'

  return
end
subroutine smb_factor ( n, xadj, iadj, perm, invprm, ixlnz, nofnz, xnzsub, &
  nzsub, maxsub, rchlnk, mrglnk )
!
!*******************************************************************************
!
!! SMB_FACTOR performs symbolic factorization on a permuted linear system.
!
!
!  Discussion:
!
!    The routine also sets up the compressed data structure for the system.
!
!  Modified:
!
!    30 April 2000
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer PERM(N), the permutation vector.
!
!    Input, integer INVPRM(N), the inverse permutation vector.
!
!  updated parameters-
!  maxsub-size of the subscript array nzsub.  on return,
!  it contains the number of subscripts used
!
!  output parameters-
!  ixlnz-index into the nonzero storage vector lnz.
!  (xnzsub,nzsub)-the compressed subscript vectors.
!  nofnz-the number of nonzeros found.
!
!  working parameters-
!  mrglnk-a vector of size N.  at the kth step,
!  mrglnk(k),mrglnk(mrglnk(k)) ,
!  is a list containing all those columns l(*,j)
!  with j less than k,such that its first off-
!  diagonal nonzero is l(k,j).  thus,the
!  nonzero structure of column l(*,k) can be found
!  by merging that of such columns l(*,j) with
!  the structure of a(*,k).
!  rchlnk-a vector of size N.  it is used to accumulate
!  the structure of each column l(*,k).  at the
!end of the kth step,
!  rchlnk(k),rchlnk(rchlnk(k)),
!  is the list of positions of nonzeros in column k
!  of the factor l.
!
  integer n
!
  integer iadj(*)
  integer i
  integer inz
  integer invprm(n)
  integer ixlnz(n+1)
  integer j
  integer jstop
  integer jstrt
  integer k
  integer knz
  integer kxsub
  integer lmax
  integer m
  integer marker(n)
  integer maxsub
  integer mrgk
  integer mrglnk(n)
  integer mrkflg
  integer nabor
  integer node
  integer nofnz
  integer np1
  integer nzbeg
  integer nzend
  integer nzsub(*)
  integer perm(n)
  integer rchlnk(n)
  integer rchm
  integer xadj(n+1)
  integer xnzsub(*)
!
!  Initialization.
!
  nzbeg = 1
  nzend = 0
  ixlnz(1) = 1

  do k = 1, n
    mrglnk(k) = 0
  end do
!
!  MARKER is used to test if mass symbolic elimination can be performed.  
!  That is, it is used to check whether the structure of the current 
!  column K being processed is completely determined by the single
!  column MRGLNK(K).
!
  marker = 0
!
!  For each column KNZ counts the number
!  of nonzeros in column K accumulated in rchlnk.
!
  np1 = n + 1

  do k = 1, n

    knz = 0
    mrgk = mrglnk(k)
    mrkflg = 0
    marker(k) = k
    if ( mrgk /= 0 ) then
      marker(k) = marker(mrgk)
    end if
    xnzsub(k) = nzend
    node = perm(k)
    jstrt = xadj(node)
    jstop = xadj(node+1)-1
    if ( jstrt > jstop ) then
      go to 160
    end if
!
!  Use RCHLNK to link through the structure of A(*,K) below diagonal.
!
    rchlnk(k) = np1

    do j = jstrt, jstop

      nabor = iadj(j)
      nabor = invprm(nabor)

      if ( nabor > k ) then

        rchm = k

   20       continue

        m = rchm
        rchm = rchlnk(m)
        if ( rchm <= nabor ) then
          go to 20
        end if

        knz = knz+1
        rchlnk(m) = nabor
        rchlnk(nabor) = rchm

        if ( marker(nabor) /= marker(k) ) then
          mrkflg = 1
        end if

      end if

    end do
!
!  Test for mass symbolic elimination
!
    lmax = 0

    if ( mrkflg /= 0 .or. mrgk == 0 ) then
      go to 40
    end if

    if ( mrglnk(mrgk) /= 0 ) then
      go to 40
    end if

    xnzsub(k) = xnzsub(mrgk) + 1
    knz = ixlnz(mrgk+1) - (ixlnz(mrgk)+1)
    go to 150
!
!  Link through each column I that affects l(*,k).
!
   40   continue

    i = k

   50   continue

    i = mrglnk(i)
    if ( i == 0 ) then
      go to 90
    end if

    inz = ixlnz(i+1)-(ixlnz(i)+1)
    jstrt = xnzsub(i)+1
    jstop = xnzsub(i)+inz

    if ( inz > lmax ) then
      lmax = inz
      xnzsub(k) = jstrt
    end if
!
!  Merge structure of l(*,i) in nzsub into RCHLNK.
!
    rchm = k

    do j = jstrt, jstop

      nabor = nzsub(j)

   70     continue

      m = rchm
      rchm = rchlnk(m)
      if ( rchm < nabor ) then
        go to 70
      end if

      if ( rchm /= nabor ) then
        knz = knz+1
        rchlnk(m) = nabor
        rchlnk(nabor) = rchm
        rchm = nabor
      end if

    end do

    go to 50
!
!  Check if subscripts duplicate those of another column...
!
   90   continue

    if ( knz == lmax ) then
      go to 150
    end if
!
!  ...or if tail of column K-1 matches head of column K.
!
    if ( nzbeg > nzend ) then
      go to 130
    end if

    i = rchlnk(k)

    do jstrt = nzbeg, nzend

      if ( nzsub(jstrt) == i ) then
        go to 110
      end if

      if ( nzsub(jstrt) > i ) then
        go to 130
      end if

    end do

    go to 130

  110   continue

    xnzsub(k) = jstrt

    do j = jstrt, nzend
      if ( nzsub(j) /= i ) then
        go to 130
      end if
      i = rchlnk(i)
      if ( i > n ) then
        go to 150
      end if
    end do

    nzend = jstrt-1
!
!  Copy the structure of l(*,k) from rchlnk
!  to the data structure (xnzsub,nzsub).
!
  130   continue

    nzbeg = nzend+1
    nzend = nzend+knz

    if ( nzend > maxsub ) then
      write ( *, * ) ' '
      write ( *, * ) 'SMB_FACTOR - Fatal error!'
      write ( *, * ) '  Insufficient storage for nonzero entries.'
      write ( *, * ) '  MAXSUB  = ', maxsub
      write ( *, * ) '  Exceeded by NZEND = ', nzend
      stop
    end if

    i = k

    do j = nzbeg, nzend
      i = rchlnk(i)
      nzsub(j) = i
      marker(i) = k
    end do

    xnzsub(k) = nzbeg
    marker(k) = k
!
!  Update the vector mrglnk.  Note column L(*,K) just found
!  is required to determine column L(*,J), where
!  L(J,K) is the first nonzero in L(*,K) below diagonal.
!
  150   continue

    if ( knz > 1 ) then
      kxsub = xnzsub(k)
      i = nzsub(kxsub)
      mrglnk(k) = mrglnk(i)
      mrglnk(i) = k
    end if

  160   continue

    ixlnz(k+1) = ixlnz(k) + knz

  end do

  nofnz = ixlnz(n) - 1
  maxsub = xnzsub(n)
  xnzsub(n+1) = xnzsub(n)

  return
end
subroutine sorts1 ( n, array )
!
!*******************************************************************************
!
!! SORTS1 ascending sorts integers using linear insertion.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  input parameter-
!  n-the size of integer array.
!
!  updated parameter-
!  array-the integer vector, which on output will be
!  in increasing order.
!
  integer n
!
  integer array(n)
  integer k
  integer l
  integer node
!
  do k = 2, n

    node = array(k)

    do l = k-1, 1, -1
      if ( array(l) <= node ) then
        go to 10
      end if
      array(l+1) = array(l)
    end do

10      continue

    array(l+1) = node

  end do

  return
end
subroutine subrcm ( xadj, iadj, mask, nsubg, subg, perm, xls, n )
!
!*******************************************************************************
!
!! SUBRCM finds the reverse cuthill mckee ordering for a given subgraph.
!
!
!  Discussion:
!
!    The subgraph may be disconnected.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer XADJ(N+1).  Information about row I is stored
!    in entries XADJ(I) through XADJ(I+1)-1 of IADJ.
!
!    Input, integer IADJ(*), the adjacency structure. 
!    For each row, it contains the column indices of the nonzero entries.
!
!  input parameters-
!
!  (nsubg,subg)-the given subgraph.  nsubg is the
!  the size of the subgraph, and subg contains
!  the nodes in it.
!
!  output parameter-
!  perm-the permutation vector. it is also used
!  temporarily to store a level structure.
!
!  working parameters-
!  mask-mask vector with all zeros.  it is used to
!  specify nodes in the subgraph.
!  xls-index to a level structure.  note that the level
!  structure is stored in part of perm.
!
  integer n
!
  integer iadj(*)
  integer i
  integer iccsze
  integer mask(n)
  integer nlvl
  integer node
  integer nsubg
  integer num
  integer perm(n)
  integer subg(n)
  integer xadj(n+1)
  integer xls(n+1)
!
  do i = 1, nsubg
    node = subg(i)
    mask(node) = 1
  end do
!
!  For each connected component in the subgraph, call FNROOT and RCM
!  for the ordering.
!
  num = 0

  do i = 1, nsubg

    node = subg(i)

    if ( mask(node) > 0 ) then

      call fnroot ( node, xadj, iadj, mask, nlvl, xls, perm(num+1), n )

      call rcm ( node, xadj, iadj, mask, perm(num+1), iccsze, xls, n )

      num = num + iccsze

      if ( num >= nsubg ) then
        return
      end if

    end if

  end do

  return
end
subroutine ts_factor ( nblks, xblk, father, diag, xenv, env, xnonz, nonz, &
  nzsubs, first, n, ierror )
!
!*******************************************************************************
!
!! TS_FACTOR performs the symmetric factorization of a tree-partitioned system.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  input parameters-
!
!  (nblks,xblk,father)-the tree partitioning.
!  xenv-the envelope index vector.
!  (xnonz,nonz,nzsubs)-the off-diagonal nonzeros in
!  the original matrix.
!
!  updated parameters-
!  (diag,env)-storage arrays for the envelope of
!  the diagonal blocks of the matrix. on output,
!  contains the diagonal blocks of the factor.
!
!  working parameter-
!  first-temporary vector used to facilitate the
!  indexing to the vector nonz (or nzsubs)
!  for non-null subcolumns in off-diagonal
!  blocks.
!
!    Output, integer IERROR, error flag.
!    0, no error, the factorization was carried out.
!    1, the matrix is not positive definite.
!
  integer nblks
  integer n
!
  integer blksze
  integer col
  integer col1
  integer colbeg
  integer colend
  integer colsze
  real diag(n)
  real env(*)
  integer father(n)
  integer first(n)
  integer fnz
  integer fnz1
  integer i
  integer ierror
  integer istop
  integer istrt
  integer isub
  integer j
  integer jstop
  integer jstrt
  integer k
  integer kenv
  integer kenv0
  integer kfathr
  real nonz(*)
  integer nzsubs(*)
  integer row
  integer rowbeg
  integer rowend
  real s
  real temp(n)
  integer xblk(nblks+1)
  integer xenv(n+1)
  integer xnonz(n+1)
!
  ierror = 0

  do i = 1, n
    temp(i) = 0.0
    first(i) = xnonz(i)
  end do
!
!  Loop through the blocks.
!
  do k = 1, nblks

    rowbeg = xblk(k)
    rowend = xblk(k+1) - 1
    blksze = rowend - rowbeg + 1

    call es_factor ( blksze, xenv(rowbeg), env, diag(rowbeg), ierror )

    if ( ierror .ne. 0 ) then
      return
    end if
!
!  Modify the father diagonal block A(FATHER(K),FATHER(K)) from the 
!  off-diagonal block A(K,FATHER(K)).
!
    kfathr = father(k)

    if ( kfathr <= 0 ) then
      go to 160
    end if

    colbeg = xblk(kfathr)
    colend = xblk(kfathr+1) - 1
!
!  Find the first and last non-null column in the off-diagonal block.
!  Reset COLBEG and COLEND.
!
    do col = colbeg, colend
      jstrt = first(col)
      jstop = xnonz(col+1)-1
      if ( jstop >= jstrt .and. nzsubs(jstrt) <= rowend ) then
        go to 30
      end if
    end do

   30   continue

    colbeg = col
    col = colend

    do col1 = colbeg, colend
      jstrt = first(col)
      jstop = xnonz(col+1)-1
      if ( jstop >= jstrt .and. nzsubs(jstrt) <= rowend ) then
        go to 50
      end if
      col = col-1
    end do

   50   continue

    colend = col

    do col = colbeg, colend

      jstrt = first(col)
      jstop = xnonz(col+1)-1
!
!  Test for null subcolumn.  FNZ stores the first nonzero subscript
!  in the block column.
!
      if ( jstop < jstrt ) then
        go to 130
      end if

      fnz = nzsubs(jstrt)
      if ( fnz > rowend ) then
        go to 130
      end if
!
!  Unpack a column in the off-diagonal block and perform upper and
!  lower solves on the unpacked column.
!
      do j = jstrt, jstop
        row = nzsubs(j)
        if ( row > rowend ) then
          go to 70
        end if
        temp(row) = nonz(j)
      end do

70        continue

      colsze = rowend - fnz + 1

      call el_solve ( colsze, xenv(fnz), env, diag(fnz), temp(fnz) )

      call eu_solve ( colsze, xenv(fnz), env, diag(fnz), temp(fnz) )
!
!  Do the modification by looping through
!  the columns and forming inner products.
!
      kenv0 = xenv(col+1)-col

      do col1 = colbeg, colend

        istrt = first(col1)
        istop = xnonz(col1+1)-1
!
!  Check to see if subcolumn is null.
!
        fnz1 = nzsubs(istrt)
        if ( istop < istrt .or. fnz1 > rowend ) then
          go to 110
        end if
!
!  Check if inner product should be done.
!
        if ( fnz1 < fnz ) then
          go to 110
        end if

        if ( fnz1 == fnz .and. col1 < col ) then
          go to 110
        end if

        s = 0.0
        do i = istrt, istop
          isub = nzsubs(i)
          if ( isub > rowend ) then
            go to 90
          end if
          s = s+temp(isub)*nonz(i)
        end do
!
!  Modify ENV or DIAG.
!
   90       continue

        if ( col1 /= col ) then

          kenv = kenv0+col1
          if ( col1 > col ) then
            kenv = xenv(col1+1)-col1+col
          end if
          env(kenv) = env(kenv)-s

        else

          diag(col1) = diag(col1)-s

        end if

  110       continue

      end do
!
!  Reset part of the TEMP vector to zero.
!
      do row = fnz, rowend
        temp(row) = 0.0
      end do

  130     continue

    end do
!
!  Update the first vector for columns in FATHER(K) block, so that it
!  will index to the beginning of the next off-diagonal block to be
!  considered.
!
    do col = colbeg, colend

      jstrt = first(col)
      jstop = xnonz(col+1)-1

      if ( jstop >= jstrt ) then

        do j = jstrt, jstop

          row = nzsubs(j)

          if ( row > rowend ) then
            first(col) = j
            go to 150
          end if

        end do

      first(col) = jstop+1

      end if

  150     continue

    end do

  160   continue

  end do

  return
end
subroutine ts_solve ( nblks, xblk, diag, xenv, env, xnonz, nonz, nzsubs, rhs, &
  n )
!
!*******************************************************************************
!
!! TS_SOLVE solves a tree-partitioned factored system by implicit back substitution.
!
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  input parameters-
!  (nblks,xblk)-the partitioning.
!  (xenv,env)-envelope of the diagonal blocks.
!  (xnonz,nonz,nzsubs)-data structure for the off-
!  block diagonal nonzeros.
!
!  updated parameters-
!  rhs-on input it contains the right hand vector.
!  on output, the solution vector.
!
  integer nblks
  integer n
!
  integer col
  integer col1
  integer col2
  real diag(n)
  real env(*)
  integer i
  integer j
  integer jstop
  integer jstrt
  integer last
  integer ncol
  real nonz(*)
  integer nrow
  integer nzsubs(*)
  real rhs(n)
  integer row
  integer row1
  integer row2
  real s
  real temp(n)
  integer xblk(nblks+1)
  integer xenv(n+1)
  integer xnonz(n+1)
!
!  Forward substitution.
!
  do i = 1, nblks

    row1 = xblk(i)
    row2 = xblk(i+1)-1
    last = xnonz(row2+1)
!
!  Modify the right hand side vector by the product of the
!  off-diagonal block with the corresponding part of the right hand side.
!
    if ( i /= 1 .and. last /= xnonz(row1) ) then

      do row = row1, row2

        jstrt = xnonz(row)

        if ( jstrt /= last ) then

          jstop = xnonz(row+1)-1

          if ( jstop >= jstrt ) then

            s = 0.0
            do j = jstrt, jstop
              col = nzsubs(j)
              s = s+rhs(col)*nonz(j)
            end do

            rhs(row) = rhs(row)-s

          end if

        end if

      end do

    end if

    nrow = row2 - row1 + 1

    call el_solve ( nrow, xenv(row1), env, diag(row1), rhs(row1) )

    call eu_solve ( nrow, xenv(row1), env, diag(row1), rhs(row1) )

  end do
!
!  Backward solution.
!
  if ( nblks == 1 ) then
    return
  end if

  last = xblk(nblks)-1

  temp(1:last) = 0.0

  i = nblks
  col1 = xblk(i)
  col2 = xblk(i+1)-1

   60 continue

  if ( i == 1 ) then
    return
  end if

  last = xnonz(col2+1)

  if ( last == xnonz(col1) ) then
    go to 90
  end if
!
!  Multiply the off-diagonal block by the corresponding
!  part of the solution vector and store in TEMP.
!
  do col = col1, col2

    s = rhs(col)

    if ( s /= 0.0 ) then

      jstrt = xnonz(col)
      if ( jstrt == last ) then
        go to 90
      end if

      jstop = xnonz(col+1)-1

      do j = jstrt, jstop
        row = nzsubs(j)
        temp(row) = temp(row)+s*nonz(j)
      end do

    end if

  end do

   90 continue

  i = i-1
  col1 = xblk(i)
  col2 = xblk(i+1)-1
  ncol = col2-col1+1

  call el_solve ( ncol,xenv(col1),env,diag(col1),temp(col1))

  call eu_solve ( ncol,xenv(col1),env,diag(col1),temp(col1))

  do j = col1, col2
    rhs(j) = rhs(j)-temp(j)
  end do

  go to 60
end
