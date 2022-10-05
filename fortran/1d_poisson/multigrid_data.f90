!----------------------------------------------------------------------

module multigrid_parameters

!------------------------------------------------------------

  integer, parameter :: max_level = 16, base_grid = 3,    &
                        pre_relax = 4, post_relax = 4,    &
                        n_vcycle = 20

  double precision, parameter :: xmin = 0.d0, xmax = 1.d0,      &
                                 pi = 4.d0*atan(1.d0),          &
                                 pi2 = pi**2,                   &
                                 rho_c = 1.28d-3,               &
                                 N = 1.d0, K = 1.d2,            &
                                 alpha2 = (N+1.d0)*K/(4.d0*pi), &
                                 alpha = sqrt(alpha2),          &
                                 r_s = pi*alpha

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

  double precision, dimension(:), allocatable :: exact

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

    allocate (exact(n), stat = err)
    if (err /= 0) then
      write(*,*) "error in memory allocation"
      stop
    endif

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

    deallocate (exact)

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

    double precision :: x, r_x, rs2, temp
    integer :: i, n

!--------------------------------------------------

    n = u%nx
    rs2 = r_s**2

    do i = 1, n

      u%q(i) = 0.d0
      exact(i) = 0.d0

    enddo

    i = 1
    u%src(i) = 4.d0*pi*rho_c*rs2
    exact(i) = -8.d0*pi*rho_c*alpha2
    do i = 2, n-1

      x = u%x(i)
      r_x = u%r_x(i)
      temp = x/(1.d0-x)
      if (x > 5.d-1) then
        u%src(i) = 0.d0
        exact(i) = -4.d0*pi*alpha2*rho_c/temp
      else
        u%src(i) = 4.d0*pi*rho_c*rs2*(sin(pi*temp)/(pi*temp))/(1.d0-x)**4
        exact(i) = -4.d0*pi*rho_c*alpha2*(1.d0+sin(pi*temp)/(pi*temp))
      endif

    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine initialize_finest_grid

!------------------------------------------------------------

end module multigrid_data

!----------------------------------------------------------------------
