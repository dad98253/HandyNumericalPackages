!  femprb.f90  02 June 2000
!
program femprb
!
!*******************************************************************************
!
!! FEMPRB calls the various FEMPACK tests.
!
  write ( *, * ) ' '
  write ( *, * ) 'FEMPRB'
  write ( *, * ) '  Begin FEMPACK tests.'

  call test01
  call test02
  call test03
  call test04

  write ( *, * ) ' '
  write ( *, * ) 'FEMPRB'
  write ( *, * ) '  Normal end of FEMPACK tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests the shape routines.
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  SHAPE_TEST tests the shape routines.'

  call shape_test ( 'Q4' )

  call shape_test ( 'Q8' )

  call shape_test ( 'Q9' )

  call shape_test ( 'Q12' )

  call shape_test ( 'Q16' )

  call shape_test ( 'QL' )

  call shape_test ( 'T3' )

  call shape_test ( 'T6' )

  call shape_test ( 'T10' )

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests the grid routines.
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  Test the grid routines.'

  call grid_test ( 'Q4' )

  call grid_test ( 'Q8' )

  call grid_test ( 'Q9' )

  call grid_test ( 'Q12' )

  call grid_test ( 'Q16' )

  call grid_test ( 'QL' )

  call grid_test ( 'T3' )

  call grid_test ( 'T6' )

  call grid_test ( 'T10' )

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests the map routines.
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  Test the map routines.'

  call map_test ( 'Q4' )

  call map_test ( 'Q8' )

  call map_test ( 'Q9' )

  call map_test ( 'Q12' )

  call map_test ( 'Q16' )

  call map_test ( 'QL' )

  call map_test ( 'T3' )

  call map_test ( 'T6' )

  call map_test ( 'T10' )

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests DIV_BIL.
!
  integer, parameter :: m = 21
  integer, parameter :: n = 13
  integer, parameter :: maxm = m
!
  real diff
  real div(maxm, n-1)
  real dudx
  real dudy
  real dvdx
  real dvdy
  integer i
  integer j
  real temp
  real u(maxm, n)
  real v(maxm, n)
  real vort(maxm, n-1)
  real x
  real xhi
  real xlo
  real xm
  real y
  real yhi
  real ylo
  real ym
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  DIV_BIL estimates divergence and vorticity.'
  write ( *, * ) ' '
  write ( *, * ) '  Original U, V data forms ', m, ' rows and '
  write ( *, * ) '  ', n, ' columns.'
!
!  Set limits of the data points.
!
  xlo = 0.0
  xhi = 1.0
  ylo = 0.0
  yhi = 2.0
!
!  Put dummy data into U and V at the data nodes.
!
  write ( *, * ) ' '
  write ( *, * ) '  I, J, X, Y, U(I,J), V(I,J)'
  write ( *, * ) ' '

  do i = 1, m
    y = ( real ( m - i ) * ylo + real ( i - 1 ) * yhi ) / real ( m - 1 )
    do j = 1, n
      x = ( real ( n - j ) * xlo + real ( j - 1 ) * xhi ) / real ( n - 1 )
      u(i,j) = x * y
      v(i,j) = sin ( x**2 + y**2 )
      write ( *, '(2i4,4g14.6)' ) i, j, x, y, u(i,j), v(i, j)
    end do
    write ( *, * ) ' '
  end do
!
!  Get DIV and VORT.
!
  call div_bil ( div, m, maxm, n, u, v, vort, xhi, xlo, yhi, ylo )
!
!  Compare computed and known values at the centers of the elements.
!
  write ( *, * ) ' '
  write ( *, * ) 'I, J, XM(I, j), YM(I, j)'
  write ( *, * ) ' '

  do i = 1, m-1

    ym = ( real(2*m-2*i-1) * ylo + real(2*i-1) * yhi ) / real ( 2 * m - 2 )

    do j = 1, n-1

      xm = ( real(2*n-2*j-1) * xlo + real(2*j-1) * xhi ) / real( 2 * n - 2 )

      write ( *, '(2i4,3g14.6)' ) i, j, xm, ym

    end do
    write ( *, * ) ' '
  end do

  write ( *, * ) ' '
  write ( *, * ) '  I, J, DIV(I, j), Exact Divergence, Difference'
  write ( *, * ) ' '

  do i = 1, m - 1

    ym = ( real(2*m-2*i-1) * ylo + real(2*i-1) * yhi ) / real ( 2 * m - 2 )

    do j = 1, n - 1

      xm = ( real(2*n-2*j-1) * xlo + real(2*j-1) * xhi ) / real ( 2 * n - 2 )

      dudx = ym
      dudy = xm
      dvdx = 2.0 * xm * cos ( xm**2 + ym**2 )
      dvdy = 2.0 * ym * cos ( xm**2 + ym**2 )

      temp = dudx + dvdy
      diff = div(i, j) - temp

      write ( *, '(2i4,3g14.6)' ) i, j, div(i,j), temp, diff

    end do
    write ( *, * ) ' '
  end do

  write ( *, * ) ' '
  write ( *, * ) '  I, J, VORT(I, j), Exact Vorticity, Difference'
  write ( *, * ) ' '

  do i = 1, m - 1

    y = ( real(2*m-2*i-1) * ylo + real(2*i-1) * yhi ) / real ( 2 * m - 2 )

    do j = 1, n-1

      x = ( real(2*n-2*j-1) * xlo + real(2*j-1) * xhi ) / real ( 2 * n - 2 )

      dudx = y
      dudy = x
      dvdx = 2.0 * x * cos ( x**2 + y**2 )
      dvdy = 2.0 * y * cos ( x**2 + y**2 )

      temp = dvdx - dudy
      diff = vort(i, j) - temp

      write ( *, '(2i4,3g14.6)' ) i, j, vort(i, j), temp, diff

    end do
    write ( *, * ) ' '
  end do

  return
end
