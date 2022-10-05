!----------------------------------------------------------------------

program multigrid_1d

!--------------------------------------------------

  use multigrid_parameters
  use multigrid_data
  use multigrid_functions

!--------------------------------------------------

  implicit none

  type (mg_data_type) :: mg_data(max_level)

  integer :: i, vcycle, level, j

  double precision :: rms_res, rms_exact
  real :: t1, t2

!--------------------------------------------------

  call allocate_data (mg_data)
  call grid (mg_data)
  call initialize_finest_grid (mg_data(1))

!--------------------------------------------------

  call cpu_time(t1)

  call write_zero (mg_data(1)%nx,mg_data(1)%q)

  do vcycle = 1, n_vcycle

!----------------------------------------
! start V-cycle
!----------------------------------------

    do level = 1, max_level-1

      do i = 1, pre_relax

        call relaxation (mg_data(level)%nx,     &
                         mg_data(level)%dx,     &
                         mg_data(level)%dx2,    &
                         mg_data(level)%r_dx,   &
                         mg_data(level)%r_dx2,  &
                         mg_data(level)%x,      &
                         mg_data(level)%r_x,    &
                         mg_data(level)%q,      &
                         mg_data(level)%src)
      enddo

      call residual (mg_data(level)%nx,     &
                     mg_data(level)%dx,     &
                     mg_data(level)%dx2,    &
                     mg_data(level)%r_dx,   &
                     mg_data(level)%r_dx2,  &
                     mg_data(level)%x,      &
                     mg_data(level)%r_x,    &
                     mg_data(level)%q,      &
                     mg_data(level)%src,    &
                     mg_data(level)%res)

      call restriction (mg_data(level+1)%nx,    &
                        mg_data(level)%nx,      &
                        mg_data(level+1)%src,   &
                        mg_data(level)%res)


      call write_zero (mg_data(level+1)%nx,   &
                       mg_data(level+1)%q)

    enddo


!----------------------------------------

    level = max_level
    call relaxation_maxlevel (mg_data(level)%nx,      &
                              mg_data(level)%dx,      &
                              mg_data(level)%dx2,     &
                              mg_data(level)%r_dx,    &
                              mg_data(level)%r_dx2,   &
                              mg_data(level)%x,       &
                              mg_data(level)%r_x,     &
                              mg_data(level)%q,       &
                              mg_data(level)%src)

!----------------------------------------

    do level = max_level-1, 1, -1

      call prolongation (mg_data(level+1)%nx,   &
                         mg_data(level)%nx,     &
                         mg_data(level+1)%q,    &
                         mg_data(level)%err)

      call add_array (mg_data(level)%nx,    &
                      mg_data(level)%q,     &
                      mg_data(level)%err)

      do i = 1, post_relax

        call relaxation (mg_data(level)%nx,       &
                         mg_data(level)%dx,       &
                         mg_data(level)%dx2,      &
                         mg_data(level)%r_dx,     &
                         mg_data(level)%r_dx2,    &
                         mg_data(level)%x,        &
                         mg_data(level)%r_x,      &
                         mg_data(level)%q,        &
                         mg_data(level)%src)

      enddo

    enddo

    level = 1

    call rms_error (mg_data(level)%nx,      &
                    mg_data(level)%dx,      &
                    mg_data(level)%dx2,     &
                    mg_data(level)%r_dx,    &
                    mg_data(level)%r_dx2,   &
                    mg_data(level)%x,       &
                    mg_data(level)%r_x,     &
                    mg_data(level)%q,       &
                    mg_data(level)%src,     &
                    exact,                  &
                    rms_res, rms_exact)

    write(*,*) vcycle, rms_res, rms_exact

!----------------------------------------
! end of V-cycle
!----------------------------------------

  enddo

  call cpu_time(t2)

  write(*,*) (t2-t1)/n_vcycle
!--------------------------------------------------

!  do i = 1, mg_data(1)%nx
!
!    write(*,*) mg_data(1)%x(i),mg_data(1)%q(i),exact(i)
!
!  enddo

!----------------------------------------

  call deallocate_data (mg_data)

!----------------------------------------------------------------------

end program multigrid_1d

!----------------------------------------------------------------------
