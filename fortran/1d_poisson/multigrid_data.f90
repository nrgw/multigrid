!----------------------------------------------------------------------

module multigrid_parameters

!------------------------------------------------------------

  integer, parameter :: max_level = 10, base_grid = 3,    &
                        pre_relax = 4, post_relax = 4,    &
                        n_vcycle = 20

  double precision, parameter :: xmin = 0.d0, xmax = 1.d0,    &
                                 pi = 4.d0*atan(1.d0),        &
                                 pi2 = pi**2,                 &
                                 rho_c = 1.28d-3,             &
                                 r_s = 8.d0

!------------------------------------------------------------

end module multigrid_parameters

!----------------------------------------------------------------------

module multigrid_data

!------------------------------------------------------------

  use multigrid_parameters

!------------------------------------------------------------

  implicit none

  type :: mg_data_type

    integer :: nx
    double precision :: dx, dx2, r_dx, r_dx2
    double precision, dimension(:), allocatable :: x, r_x,            &
                                                   q, src, res, err

  end type mg_data_type

!------------------------------------------------------------

  contains

!------------------------------------------------------------

  subroutine allocate_data (u)

!--------------------------------------------------

    implicit none

    type (mg_data_type) :: u(max_level)

    integer :: l, n, err

!--------------------------------------------------

    n = base_grid
    u(max_level)%nx = n
    allocate (u(max_level)%x(n),      &
              u(max_level)%r_x(n),    &
              u(max_level)%q(n),      &
              u(max_level)%src(n),    &
              u(max_level)%res(n),    &
              u(max_level)%err(n),    &
              stat = err)

    if (err /= 0) then
      write(*,*) "error in memory allocation"
      stop
    endif

    do l = max_level-1, 1, -1

      n = (n-1)*2+1
      u(l)%nx = n
      allocate (u(l)%x(n),      &
                u(l)%r_x(n),    &
                u(l)%q(n),      &
                u(l)%src(n),    &
                u(l)%res(n),    &
                u(l)%err(n),    &
                stat = err)

      if (err /= 0) then
        write(*,*) "error in memory allocation"
        stop
      endif


    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine allocate_data

!------------------------------------------------------------

  subroutine deallocate_data (u)

!--------------------------------------------------

    implicit none

    type (mg_data_type) :: u(max_level)

    integer :: l

!--------------------------------------------------

    do l = 1, max_level

      deallocate (u(l)%x,     &
                  u(l)%r_x,   &
                  u(l)%q,     &
                  u(l)%src,   &
                  u(l)%res,   &
                  u(l)%err)

    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine deallocate_data

!------------------------------------------------------------

  subroutine grid (u)

!--------------------------------------------------

    implicit none

    type (mg_data_type) :: u(max_level)

    double precision :: dx, dx2, r_dx, r_dx2
    integer :: l, n, i

!--------------------------------------------------

    do l = 1, max_level

      n = u(l)%nx
      dx = (xmax-xmin)/dble(n-1)
      dx2 = dx**2
      r_dx = 1.d0/dx
      r_dx2 = r_dx**2

      u(l)%dx = dx
      u(l)%dx2 = dx2
      u(l)%r_dx = r_dx
      u(l)%r_dx2 = r_dx2

      do i = 1, n

        u(l)%x(i) = xmin+dx*(i-1)

      enddo


      u(l)%r_x(1) = 0.d0
      do i = 2, n

        u(l)%r_x(i) = 1.d0/u(l)%x(i)

      enddo

     

    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine grid

!------------------------------------------------------------

  subroutine initialize_finest_grid (u)

!--------------------------------------------------

    implicit none

    type (mg_data_type) :: u

    double precision :: x, rs2
    integer :: i, n

!--------------------------------------------------

    n = u%nx
    rs2 = r_s**2

    do i = 1, n

      u%q(i) = 0.d0

    enddo

    do i = 1, n-1

      x = u%x(i)
      if (x > 5.d-1) then
        u%src(i) = 0.d0
      else
        u%src(i) = 4.d0*pi*rho_c*rs2*(1.d0-(x/(1.d0-x))**2)/(1.d0-x)**4
      endif

    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine initialize_finest_grid

!------------------------------------------------------------

end module multigrid_data

!----------------------------------------------------------------------
