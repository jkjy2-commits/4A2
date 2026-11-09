      
      subroutine set_secondary(av,g)

!     This subroutine calculates the secondary flow variables from the primary
!     ones at every node

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g

!     Define any further variables you may need
!     INSERT
      real, dimension(g%ni,g%nj) :: t

!     The primary flow variables are "ro", "roe", "rovx" and "rovy", these are 
!     the conserved quantities within the Euler equations. Write the code to
!     calculate the secondary flow variables which are the velocity components
!     "vx" and "vy", the static pressure "p" and the stagnation enthalpy
!     "hstag". These are needed at every timestep, there is no need for any 
!     loops as the operations can be performed elementwise, although you may
!     wish to define some intermediate variables to improve readability.
!     INSERT
      g%vx = g%rovx / g%ro
      g%vy = g%rovy / g%ro

      t = (g%roe / g%ro - 0.5 * (g%vx**2 + g%vy**2))/ av%cv

      g%p = g%ro * av%rgas * t

      g%hstag = av%cp * t + &
            0.5 * (g%vx**2 + g%vy**2)

      end subroutine set_secondary


