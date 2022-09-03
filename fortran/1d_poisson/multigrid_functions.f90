!----------------------------------------------------------------------

module multigrid_functions

!------------------------------------------------------------

  implicit none

  double precision, parameter :: one_over_three = 1.d0/3.d0,          &
                                 one_over_six = 5.d-1*one_over_three

!------------------------------------------------------------

  contains

!------------------------------------------------------------

  subroutine relaxation (nx, dx, dx2, r_dx, r_dx2, x, r_x, q, s)

!--------------------------------------------------

    implicit none
    
    integer :: nx
    double precision :: dx, dx2, r_dx, r_dx2
    double precision :: x(nx), r_x(nx), q(nx), s(nx)
    
    integer :: i, j

!--------------------------------------------------
! Gauss-Seidel red black relaxation
!--------------------------------------------------

    i = 1
    q(i) = q(i+1)-dx2*one_over_six*s(i)
    do j = 3, 2, -1
      do i = j, nx-1, 2
        
        q(i) = 5.d-1*(q(i+1)+q(i-1)   &
                     +dx*r_x(i)*(q(i+1)-q(i-1))    &
                     -dx2*s(i))

      enddo
    enddo

!--------------------------------------------------

  end subroutine relaxation

!------------------------------------------------------------

  subroutine residual (nx, dx, dx2, r_dx, r_dx2, x, r_x, q, s, r)

!--------------------------------------------------

    implicit none
    
    integer :: nx
    double precision :: dx, dx2, r_dx, r_dx2
    double precision :: x(nx), r_x(nx), q(nx), s(nx), r(nx)
    
    integer :: i

!--------------------------------------------------

    r(nx) = 0.d0

    i = 1
    r(i) = s(i)-6.d0*(q(i+1)-q(i))*r_dx2
    do i = 2, nx-1

      r(i) = s(i)-(q(i+1)-2.d0*q(i)+q(i-1))*r_dx2   &
            -r_x(i)*(q(i+1)-q(i-1))*r_dx
    
    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine residual

!------------------------------------------------------------

  subroutine rms_error (nx, dx, dx2, r_dx, r_dx2, x, r_x, q, s, exact, &
                        rms_res, rms_exact)

!--------------------------------------------------

    implicit none
    
    integer :: nx
    double precision :: dx, dx2, r_dx, r_dx2
    double precision :: x(nx), r_x(nx), q(nx), s(nx), exact(nx)
    double precision :: rms_res, rms_exact

    double precision :: temp1, temp2
    integer :: i

!--------------------------------------------------

    rms_exact = 0.d0

    i = 1
    temp1 = s(i)-6.d0*(q(i+1)-q(i))*r_dx2
    rms_res = temp1**2

    do i = 2, nx-1

      temp1 = (s(i)-(q(i+1)-2.d0*q(i)+q(i-1))*r_dx2   &
              -r_x(i)*(q(i+1)-q(i-1))*r_dx)

      rms_res = rms_res+temp1**2

    enddo

    do i = 1, nx
      rms_exact = rms_exact+(q(i)-exact(i))**2
    enddo
    
    temp2 = 1.d0/dble(nx)
    rms_res = sqrt(rms_res*temp2)
    rms_exact = sqrt(rms_exact*temp2)

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine rms_error

!------------------------------------------------------------

  subroutine relaxation_maxlevel (nx, dx, dx2, r_dx, r_dx2,   &
                                  x, r_x, q, s)

!--------------------------------------------------

    implicit none
    
    integer :: nx
    double precision :: dx, dx2, r_dx, r_dx2
    double precision :: x(nx), r_x(nx), q(nx), s(nx)

    integer :: i

!--------------------------------------------------

    q(1) = q(3)-dx2/(1.d0+dx*r_x(2))*(one_over_three*s(1)+s(2))
    q(2) = q(1)+one_over_six*dx2*s(1)

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine relaxation_maxlevel

!------------------------------------------------------------

  subroutine restriction (n_coarse, n_fine, q_coarse, q_fine)

!--------------------------------------------------

    implicit none

    integer :: n_coarse, n_fine
    double precision :: q_fine(n_fine)
    double precision :: q_coarse(n_coarse)

    integer :: i, i_fine

!--------------------------------------------------

    q_coarse(n_coarse) = 0.d0

    q_coarse(1) = 5.d-1*(q_fine(1)+q_fine(2)) !reflective bc
    do i = 2, n_coarse-1

      i_fine = 2*i-1
      q_coarse(i) = 2.5d-1*(q_fine(i_fine-1)     &
                           +q_fine(i_fine+1)     &
                           +2.d0*q_fine(i_fine))

    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine restriction

!------------------------------------------------------------

  subroutine prolongation (n_coarse,n_fine,q_coarse,q_fine)

!--------------------------------------------------

    implicit none

    integer :: n_coarse, n_fine
    double precision :: q_coarse(n_coarse)
    double precision :: q_fine(n_fine)

    integer :: i, i_fine

!--------------------------------------------------

    do i = 1, n_coarse

      i_fine = 2*i-1
      q_fine(i_fine) = q_coarse(i)
    
    enddo
    
    do i = 2, n_fine, 2
      
      q_fine(i) = 5.d-1*(q_fine(i-1)+q_fine(i+1))
    
    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine prolongation

!------------------------------------------------------------

  subroutine write_zero (n,q)

!--------------------------------------------------

    implicit none
    
    integer :: n
    double precision :: q(n)

    integer :: i

!--------------------------------------------------

    do i = 1, n
      
      q(i) = 0.d0
    
    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine write_zero

!------------------------------------------------------------

  subroutine copy_array (n,from,to)

!--------------------------------------------------

    implicit none
    
    integer :: n
    double precision :: from(n)
    double precision :: to(n)

    integer :: i

!--------------------------------------------------

    do i = 1, n
      
      to(i) = from(i)
    
    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine copy_array

!------------------------------------------------------------

  subroutine add_array (n,a,b)

!--------------------------------------------------

    implicit none
    
    integer :: n
    double precision :: a(n)
    double precision :: b(n)
    
    integer :: i

!--------------------------------------------------

    do i = 1, n
      
      a(i) = a(i)+b(i)
    
    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine add_array

!------------------------------------------------------------

  subroutine subtract_array (n,a,b)

!--------------------------------------------------

    implicit none
    
    integer :: n
    double precision :: a(n)
    double precision :: b(n)
    
    integer :: i

!--------------------------------------------------

    do i = 1, n
      
      a(i) = a(i)-b(i)
    
    enddo

!--------------------------------------------------

    return

!--------------------------------------------------

  end subroutine subtract_array

!------------------------------------------------------------

end module multigrid_functions

!----------------------------------------------------------------------

