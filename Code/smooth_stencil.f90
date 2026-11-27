
      module smooth_stencil

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_array(av,prop,prop_ref)

!     This subroutine smooths "prop" to stabilise the calculation, the basic 
!     solver uses second order smoothing, many improvements are possible.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(inout) :: prop(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg_2, prop_avg_4
      real, intent(in) :: prop_ref
      real, dimension(size(prop,1),size(prop,2)) :: avg_4_i, avg_4_j, sfac_loc
      real, dimension(size(prop,1),size(prop,2)) :: sf2, sf4
      integer :: ni, nj

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)
!     IMPROVE


!     Calculate the average values at the nodes in the interior region of the
!     mesh, use the four neighbouring nodes in the plus and minus i and 
!     j-directions.
!     INSERT
      prop_avg_2(2:ni-1,2:nj-1) = 0.25 * (prop(1:ni-2,2:nj-1) + prop(3:ni,2:nj-1) + &
          prop(2:ni-1,1:nj-2) + prop(2:ni-1,3:nj))

!     Edge values are also averaged in both the i and j-directions. Parallel to
!     the boundary the averaging is centred, the averages of two nodes are taken
!     either side of the current point. Perpendicular to the boundary the
!     algorithm is one-sided, the value at the current point is extrapolated
!     from the values at two nodes away from the boundary point.
!     INSERT
      prop_avg_2(1,2:nj-1) = (1.0/3.0) * (prop(1,1:nj-2) + prop(1,3:nj) + &
          2* prop(2,2:nj-1) - prop(3,2:nj-1))
      prop_avg_2(ni,2:nj-1) = (1.0/3.0) * (prop(ni,1:nj-2) + prop(ni,3:nj) + &
          2* prop(ni-1,2:nj-1) - prop(ni-2,2:nj-1))
      prop_avg_2(2:ni-1,1) = (1.0/3.0) * (prop(1:ni-2,1) + prop(3:ni,1) + &
          2* prop(2:ni-1,2) - prop(2:ni-1,3))
      prop_avg_2(2:ni-1,nj) = (1.0/3.0) * (prop(1:ni-2,nj) + prop(3:ni,nj) + &
          2* prop(2:ni-1,nj-1) - prop(2:ni-1,nj-2))   

!     The corner values are not currently smoothed
!     IMPROVE
      prop_avg_2([1,ni],[1,nj]) = prop([1,ni],[1,nj])

!     fourth order smoothing
      prop_avg_4 = prop_avg_2
!     IMPROVE
!     Cubic through i-2, i-1, i+1, i+2, evaluated at i
      avg_4_i(3:ni-2,3:nj-2) = (2.0/3.0) * (prop(2:ni-3,3:nj-2) + prop(4:ni-1,3:nj-2)) &
                            - (1.0/6.0) * (prop(1:ni-4,3:nj-2) + prop(5:ni,3:nj-2))

!     Cubic through j-2, j-1, j+1, j+2, evaluated at j
      avg_4_j(3:ni-2,3:nj-2) = (2.0/3.0) * (prop(3:ni-2,2:nj-3) + prop(3:ni-2,4:nj-1)) &
                            - (1.0/6.0) * (prop(3:ni-2,1:nj-4) + prop(3:ni-2,5:nj))

!     Combine the two directions
      prop_avg_4(3:ni-2,3:nj-2) = 0.5 * (avg_4_i(3:ni-2,3:nj-2) + avg_4_j(3:ni-2,3:nj-2))

!     Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!     take (1-sfac) * the calculated value of the property + sfac * the average 
!     of the surrounding values. 
!     INSERT
!     IMPROVE
      sfac_loc = min(av%sfac * abs(prop - prop_avg_2) / prop_ref, 0.5)
      sf4 = sfac_loc
      sf2 = 0.25 * sfac_loc
      prop = (1.0 - sf2 - sf4) * prop + sf2 * prop_avg_2 + sf4 * prop_avg_4

      end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module smooth_stencil


